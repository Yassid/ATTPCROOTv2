#!/usr/bin/env python3
"""Animate ONE AT-TPC event through the reconstruction chain with VTK, from export_event_vtk_dp.C's JSON.

  hits (drift-time order) -> pattern recognition (hits coloured by track) -> hit clustering
  -> Kalman fit (trajectory grows from the vertex, clusters light up as they are passed) -> result

  ~/venvs/vtk/bin/python vtk_event_anim.py ev.json out.mp4 [out.gif] [--rescue]

--rescue: after clustering, show the FULL-LENGTH fit (fitFull, gf_dp_cateloss) failing, FIND cutting the
track back from the far end, and the refit on the kept prefix (fit, gf_dp_find805). The rung sequence
shown (100 % fails, 80 % fails, kept %) follows find805's ladder (full, then 80, 75, ... 5); the kept
fraction is stored fit points / clusters, which FIND does not record directly.

WHAT THE KALMAN STAGE SHOWS: the SMOOTHED genfit trajectory stored in the file, revealed point by
point from the vertex. It is not a replay of the filter's internal predict/update states (those
are not stored); the growing curve is a visual metaphor for the filter walking the measurements.
DISPLAY ONLY: the stored states are walked along the path from the vertex (stored order is the
z-sort, scrambled inside a z bucket) and joined by a spline, so gaps in the measurements curve
instead of cutting chords. No number shown is computed here; all come from the exported fit.
Frames: hits/clusters in the TPC frame; the exporter has already flipped the fit z (1000 - z_fit).
"""
import json
import sys

import imageio.v2 as imageio
import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, vtk_to_numpy

RESCUE = "--rescue" in sys.argv   # full-length fit fails -> FIND shortens -> refit (needs fitFull)
argv = [a for a in sys.argv if a != "--rescue"]
src, out_mp4 = argv[1], argv[2]
out_gif = argv[3] if len(argv) > 3 else None
D = json.load(open(src))
W, H, FPS = 1280, 720, 25

hits = np.array(D["hits"], float)
hits = hits[np.argsort(hits[:, 2])]                    # drift-time order: readout sweep
tracks = D["tracks"]
trk = next(t for t in tracks if t["id"] == D["tid"])
fit = D["fit"]
sm = np.array(fit["smoothed"], float)[:, :3]
cl = np.array(trk["clusters"], float)
vtx = np.array(fit["vertex"][:3], float)


def path_order(p, start):
    """nearest-neighbour walk from the point closest to `start`"""
    n = len(p)
    used = np.zeros(n, bool)
    cur = int(((p - start) ** 2).sum(1).argmin())
    out = [cur]
    used[cur] = True
    for _ in range(n - 1):
        d = ((p - p[cur]) ** 2).sum(1)
        d[used] = np.inf
        cur = int(d.argmin())
        used[cur] = True
        out.append(cur)
    return np.array(out)


# The smoothed states come in MEASUREMENT order (production z-sort), which is scrambled inside a
# z bucket; joined in that order the curve throws chords across the helix. For DISPLAY, walk them
# along the path from the vertex -- the trajectory itself is unchanged.
sm = sm[path_order(sm, vtx)]

if RESCUE:
    full = D["fitFull"]
    vtxF = np.array(full["vertex"][:3], float)
    smF = np.array(full["smoothed"], float)[:, :3]
    smF = smF[path_order(smF, vtxF)]
    # clusters in path order from the (good) vertex; the kept prefix = clusters the FIND fit passes
    cl = cl[path_order(cl[:, :3], vtx)]
    dk = np.sqrt(((cl[:, None, :3] - sm[None, :, :]) ** 2).sum(-1).min(1))
    kept = dk < 12.0
    n_keep = int(kept.sum())
    keep_pct = 100.0 * len(sm) / len(cl)
    cl_idx = np.arange(len(cl))
else:
    # each cluster's position along the fit = index of its nearest smoothed point
    d2 = ((cl[:, None, :3] - sm[None, :, :]) ** 2).sum(-1)
    cl_idx = d2.argmin(1)
    cl = cl[np.argsort(cl_idx)]
    cl_idx = np.sort(cl_idx)

# track membership of every hit (exact match on the exported coordinates)
key = lambda p: (round(p[0], 1), round(p[1], 1), round(p[2], 1))
member = {}
for t in tracks:
    for h in t["hits"]:
        member[key(h)] = t["id"]
hit_trk = np.array([member.get(key(h), -1) for h in hits])

# ---------------------------------------------------------------- scene
ren = vtk.vtkRenderer()
ren.SetBackground(0.06, 0.07, 0.09)
ren.SetBackground2(0.16, 0.18, 0.22)
ren.GradientBackgroundOn()
rw = vtk.vtkRenderWindow()
rw.SetOffScreenRendering(1)
rw.AddRenderer(ren)
rw.SetSize(W, H)
rw.SetMultiSamples(8)


def actor_for(poly, color=None, size=None, opacity=1.0, scalars=False, lut=None, rng=None):
    m = vtk.vtkPolyDataMapper()
    m.SetInputData(poly) if isinstance(poly, vtk.vtkPolyData) else m.SetInputConnection(poly)
    if scalars:
        m.SetScalarModeToUsePointData()
        m.SetColorModeToMapScalars()
        m.SetLookupTable(lut)
        m.SetScalarRange(*rng)
        m.ScalarVisibilityOn()
    else:
        m.ScalarVisibilityOff()
    a = vtk.vtkActor()
    a.SetMapper(m)
    p = a.GetProperty()
    if color:
        p.SetColor(*color)
    if size:
        p.SetPointSize(size)
        p.RenderPointsAsSpheresOn()
    p.SetOpacity(opacity)
    ren.AddActor(a)
    return a, m


def points_poly(xyz, scal=None):
    pts = vtk.vtkPoints()
    pts.SetData(numpy_to_vtk(np.ascontiguousarray(xyz, dtype=np.float64), deep=True))
    pd = vtk.vtkPolyData()
    pd.SetPoints(pts)
    verts = vtk.vtkCellArray()
    pd.SetVerts(verts)
    if scal is not None:
        arr = numpy_to_vtk(np.ascontiguousarray(scal, dtype=np.float64), deep=True)
        arr.SetName("s")
        pd.GetPointData().SetScalars(arr)
    return pd


def set_count(pd, n):
    """show the first n points of a vertex polydata"""
    ca = vtk.vtkCellArray()
    if n > 0:
        ids = np.arange(n, dtype=np.int64)
        off = np.arange(n + 1, dtype=np.int64)
        ca.SetData(numpy_to_vtk(off, deep=True, array_type=vtk.VTK_ID_TYPE),
                   numpy_to_vtk(ids, deep=True, array_type=vtk.VTK_ID_TYPE))
    pd.SetVerts(ca)
    pd.Modified()


def set_glyph_count(pd, xyz, n):
    """vtkGlyph3D glyphs EVERY point and ignores vertex cells, so reveal by truncating the points"""
    pts = vtk.vtkPoints()
    if n > 0:
        pts.SetData(numpy_to_vtk(np.ascontiguousarray(xyz[:n], dtype=np.float64), deep=True))
    pd.SetPoints(pts)
    pd.Modified()


# detector: drift volume outline + beam axis
cyl = vtk.vtkCylinderSource()
cyl.SetRadius(292.0)
cyl.SetHeight(1000.0)
cyl.SetResolution(96)
cyl.CappingOff()
tf = vtk.vtkTransform()
tf.Translate(0, 0, 500)
tf.RotateX(90)
tfp = vtk.vtkTransformFilter()
tfp.SetTransform(tf)
tfp.SetInputConnection(cyl.GetOutputPort())
det, _ = actor_for(tfp.GetOutputPort(), color=(0.55, 0.62, 0.75), opacity=0.06)
det.GetProperty().BackfaceCullingOn()
KEEP = []  # VTK sources referenced only through an output port are freed by Python -> segfault
ring_actors = []
for z in (0.0, 1000.0):
    rs = vtk.vtkRegularPolygonSource()
    rs.SetNumberOfSides(120)
    rs.SetRadius(292.0)
    rs.SetCenter(0, 0, z)
    rs.GeneratePolygonOff()
    KEEP.append(rs)
    a, _ = actor_for(rs.GetOutputPort(), color=(0.55, 0.62, 0.75), opacity=0.5)
    a.GetProperty().SetLineWidth(1.5)
    ring_actors.append(a)
beam = vtk.vtkLineSource()
beam.SetPoint1(0, 0, -60)
beam.SetPoint2(0, 0, 1060)
beam_a, _ = actor_for(beam.GetOutputPort(), color=(0.95, 0.85, 0.35), opacity=0.55)
beam_a.GetProperty().SetLineStipplePattern(0xF0F0)
beam_a.GetProperty().SetLineWidth(2)

# hits coloured by log charge
q = np.log10(np.clip(hits[:, 3], 1, None))
lut = vtk.vtkLookupTable()
lut.SetNumberOfTableValues(256)
lut.Build()
import colorsys  # noqa: E402  (small inline colormap, no matplotlib here)
for i in range(256):
    t = i / 255.0
    r, g, b = colorsys.hsv_to_rgb(0.66 - 0.66 * t, 0.85, 0.35 + 0.65 * t)
    lut.SetTableValue(i, r, g, b, 1.0)
hit_pd = points_poly(hits[:, :3], q)
hit_a, hit_m = actor_for(hit_pd, size=5.0, scalars=True, lut=lut, rng=(np.percentile(q, 2), np.percentile(q, 98)))
set_count(hit_pd, 0)

# hits coloured by track (pattern recognition); unassigned grey
TRK_COLORS = [(0.20, 0.80, 1.00), (1.0, 0.55, 0.15), (0.6, 1.0, 0.4), (1.0, 0.4, 0.8)]
rgb = np.zeros((len(hits), 3), np.uint8)
for i, t in enumerate(hit_trk):
    c = (0.45, 0.45, 0.45) if t < 0 else TRK_COLORS[t % len(TRK_COLORS)]
    rgb[i] = (np.array(c) * 255).astype(np.uint8)
hit2_pd = points_poly(hits[:, :3])
carr = numpy_to_vtk(rgb, deep=True)
carr.SetName("rgb")
hit2_pd.GetPointData().SetScalars(carr)
set_count(hit2_pd, len(hits))
hit2_a, hit2_m = actor_for(hit2_pd, size=4.5)
hit2_m.ScalarVisibilityOn()
hit2_m.SetColorModeToDirectScalars()
hit2_a.GetProperty().SetOpacity(0.0)

# clusters: dim (formed) and lit (measured by the filter)


def sphere_glyph(radius):
    s = vtk.vtkSphereSource()
    KEEP.append(s)
    s.SetRadius(radius)
    s.SetThetaResolution(14)
    s.SetPhiResolution(14)
    return s


cl_pd = points_poly(cl[:, :3])
set_glyph_count(cl_pd, cl[:, :3], 0)
g1 = vtk.vtkGlyph3D()
g1.SetSourceConnection(sphere_glyph(2.6).GetOutputPort())
g1.SetInputData(cl_pd)
g1.ScalingOff()
cl_a, _ = actor_for(g1.GetOutputPort(), color=(0.75, 0.80, 0.90), opacity=0.0)
cl_a.GetProperty().SetSpecular(0.4)

lit_pd = points_poly(cl[:, :3])
set_glyph_count(lit_pd, cl[:, :3], 0)
g2 = vtk.vtkGlyph3D()
g2.SetSourceConnection(sphere_glyph(3.0).GetOutputPort())
g2.SetInputData(lit_pd)
g2.ScalingOff()
lit_a, _ = actor_for(g2.GetOutputPort(), color=(1.0, 0.78, 0.25))
lit_a.GetProperty().SetSpecular(0.6)
lit_a.GetProperty().SetAmbient(0.25)

# fitted trajectory tube, grown from the vertex
traj_pts = vtk.vtkPoints()
traj_pd = vtk.vtkPolyData()
traj_pd.SetPoints(traj_pts)
# spline through the stored states: across a gap in the measurements (the track passing near the
# beam axis, where there are no clusters) a polyline cuts a straight chord; the spline curves.
spl = vtk.vtkSplineFilter()
spl.SetInputData(traj_pd)
spl.SetSubdivideToLength()
spl.SetLength(1.5)
tube = vtk.vtkTubeFilter()
tube.SetInputConnection(spl.GetOutputPort())
tube.SetRadius(1.6)
tube.SetNumberOfSides(16)
tube_a, _ = actor_for(tube.GetOutputPort(), color=(1.0, 0.25, 0.30))
tube_a.GetProperty().SetSpecular(0.5)
tube_a.VisibilityOff()


def set_traj(n):
    traj_pts.Reset()
    pl = vtk.vtkCellArray()
    if n >= 2:
        for p in sm[:n]:
            traj_pts.InsertNextPoint(*p)
        pl.InsertNextCell(n)
        for i in range(n):
            pl.InsertCellPoint(i)
    traj_pd.SetLines(pl)
    traj_pd.Modified()
    tube_a.SetVisibility(n >= 2)


def make_tube(color, curve):
    pts = vtk.vtkPoints()
    pd = vtk.vtkPolyData()
    pd.SetPoints(pts)
    sp = vtk.vtkSplineFilter()
    sp.SetInputData(pd)
    sp.SetSubdivideToLength()
    sp.SetLength(1.5)
    tb = vtk.vtkTubeFilter()
    tb.SetInputConnection(sp.GetOutputPort())
    tb.SetRadius(1.6)
    tb.SetNumberOfSides(16)
    KEEP.extend([sp, tb])
    a, _ = actor_for(tb.GetOutputPort(), color=color)
    a.GetProperty().SetSpecular(0.5)
    a.VisibilityOff()

    def grow(n):
        pts.Reset()
        pl = vtk.vtkCellArray()
        if n >= 2:
            for p in curve[:n]:
                pts.InsertNextPoint(*p)
            pl.InsertNextCell(n)
            for i in range(n):
                pl.InsertCellPoint(i)
        pd.SetLines(pl)
        pd.Modified()
        a.SetVisibility(n >= 2)
    return a, grow


if RESCUE:
    full_a, grow_full = make_tube((0.25, 0.55, 1.0), smF)
    vtxF_src = sphere_glyph(7.0)
    vtxF_src.SetCenter(*vtxF)
    vtxF_a, _ = actor_for(vtxF_src.GetOutputPort(), color=(0.25, 0.55, 1.0))
    vtxF_a.VisibilityOff()

state_src = sphere_glyph(6.5)
state_src.SetCenter(*sm[0])
state_a, _ = actor_for(state_src.GetOutputPort(), color=(1.0, 1.0, 1.0), opacity=0.9)
state_a.VisibilityOff()

vtx_src = sphere_glyph(7.0)
vtx_src.SetCenter(*vtx)
vtx_a, _ = actor_for(vtx_src.GetOutputPort(), color=(0.30, 1.0, 0.45))
vtx_a.VisibilityOff()


# text overlays
def text(x, y, size, bold=False, color=(0.92, 0.93, 0.95), just=0):
    t = vtk.vtkTextActor()
    p = t.GetTextProperty()
    p.SetFontSize(size)
    p.SetColor(*color)
    p.SetBold(bold)
    p.SetFontFamilyToArial()
    p.SetJustification(just)
    t.SetDisplayPosition(x, y)
    ren.AddViewProp(t)
    return t


title = text(28, H - 52, 30, bold=True)
title.SetInput("AT-TPC  16C(d,p)17C")
sub = text(28, H - 84, 18, color=(0.70, 0.74, 0.80))
sub.SetInput(f"run_{D['run']}  event {D['entry']}  track {D['tid']}")
stage = text(28, 34, 26, bold=True, color=(1.0, 0.82, 0.35))
detail = text(28, 10, 17, color=(0.78, 0.80, 0.85))
result = text(W - 28, H - 52, 20, color=(0.95, 0.95, 0.95), just=2)
result.GetTextProperty().SetVerticalJustificationToTop()
result.SetDisplayPosition(W - 28, H - 30)

# camera framed on the event, orbiting
lo, hi = hits[:, :3].min(0), hits[:, :3].max(0)
ctr = 0.5 * (lo + hi)
cam = ren.GetActiveCamera()
cam.SetFocalPoint(*ctr)
cam.SetPosition(ctr[0] + 900, ctr[1] + 520, ctr[2] - 300)
cam.SetViewUp(0, 1, 0)
cam.SetViewAngle(30)
ren.ResetCameraClippingRange()

light = vtk.vtkLight()
light.SetLightTypeToCameraLight()
light.SetPosition(0.4, 0.6, 1)
ren.AddLight(light)

w2i = vtk.vtkWindowToImageFilter()
w2i.SetInput(rw)
w2i.SetInputBufferTypeToRGB()
w2i.ReadFrontBufferOff()

frames = []


def grab():
    rw.Render()
    w2i.Modified()
    w2i.Update()
    im = w2i.GetOutput()
    a = vtk_to_numpy(im.GetPointData().GetScalars()).reshape(H, W, 3)[::-1]
    frames.append(a.copy())


ease = lambda t: 0.5 - 0.5 * np.cos(np.pi * np.clip(t, 0, 1))
ORBIT = 0.22  # deg / frame


def step(n, fn):
    for i in range(n):
        fn(i / max(n - 1, 1))
        cam.Azimuth(ORBIT)
        ren.ResetCameraClippingRange()
        grab()


# ---------------------------------------------------------------- timeline
nh, nc, ns = len(hits), len(cl), len(sm)
stage.SetInput("")
step(int(1.2 * FPS), lambda t: None)

stage.SetInput("1  Raw hits")
detail.SetInput(f"{nh} hits read out in drift time, coloured by charge")
step(int(4.0 * FPS), lambda t: set_count(hit_pd, int(ease(t) * nh)))
step(int(0.8 * FPS), lambda t: None)

stage.SetInput("2  Pattern recognition")
n_in = int((hit_trk >= 0).sum())
detail.SetInput(f"hits grouped into tracks: {len(tracks)} track{'s' if len(tracks) > 1 else ''}, "
                f"{n_in} of {nh} hits assigned")


def s2(t):
    hit2_a.GetProperty().SetOpacity(ease(t))
    hit_a.GetProperty().SetOpacity(1 - ease(t))


step(int(2.2 * FPS), s2)
step(int(0.6 * FPS), lambda t: None)

stage.SetInput("3  Hit clustering")
detail.SetInput(f"{n_in} hits -> {nc} charge-weighted clusters along the track")


def s3(t):
    set_glyph_count(cl_pd, cl[:, :3], int(ease(t) * nc))
    cl_a.GetProperty().SetOpacity(min(1.0, 3 * t))
    hit2_a.GetProperty().SetOpacity(1 - 0.85 * ease(t))


step(int(3.2 * FPS), s3)
step(int(0.6 * FPS), lambda t: None)

if RESCUE:
    nsF = len(smF)
    mp = 938.272088
    pS = 299.792458 * D["brho"]
    keS = np.sqrt(pS * pS + mp * mp) - mp
    rF = full.get("rvtx", np.hypot(vtxF[0], vtxF[1]))   # exporter's unrounded radius
    rK = fit.get("rvtx", np.hypot(vtx[0], vtx[1]))
    red = (1.0, 0.35, 0.35)

    stage.SetInput("4  Kalman fit on the FULL track")
    state_a.VisibilityOn()

    def r4(t):
        k = max(2, int(ease(t) * nsF))
        grow_full(k)
        state_src.SetCenter(*smF[k - 1])
        set_glyph_count(lit_pd, cl[:, :3], int(round(k / nsF * nc)))
        detail.SetInput(f"measurement {k} / {nsF}   the far end of the track pulls the helix off")

    step(int(6.0 * FPS), r4)
    state_a.VisibilityOff()
    vtxF_a.VisibilityOn()
    detail.SetInput("")
    result.GetTextProperty().SetColor(0.45, 0.68, 1.0)
    result.SetInput(f"FULL LENGTH\nKE  {full['ke']:.2f} MeV\nchi2/ndf  {full['chi2ndf']:.1f}\n"
                    f"vertex  {rF:.0f} mm off axis")
    stage.SetInput("4  Full-length fit FAILS")
    stage.GetTextProperty().SetColor(*red)
    step(int(3.0 * FPS), lambda t: None)

    stage.GetTextProperty().SetColor(1.0, 0.82, 0.35)
    stage.SetInput("5  FIND: cut the track back from the far end")
    result.SetInput("")
    vtxF_a.VisibilityOff()
    cl_a.GetProperty().SetColor(0.38, 0.40, 0.45)
    cl_a.GetProperty().SetOpacity(1.0)
    set_glyph_count(cl_pd, cl[:, :3], nc)

    def fade_full(t):
        full_a.GetProperty().SetOpacity(1 - 0.85 * ease(t))

    step(int(1.2 * FPS), fade_full)
    for label, frac, verdict in (("100 %", 1.0, "chi2/ndf %.1f  ->  fails" % full["chi2ndf"]),
                                 ("80 %", 0.80, "no pass  ->  shorter"),
                                 (f"{keep_pct:.0f} %", n_keep / nc,
                                  "chi2/ndf %.3f < 0.1  ->  longest passing prefix" % fit["chi2ndf"])):
        n_to = int(round(frac * nc))
        n_from = int(lit_pd.GetNumberOfPoints())
        detail.SetInput(f"prefix {label} of the clusters:  {verdict}")
        step(int(1.6 * FPS), lambda t, a=n_from, b=n_to: set_glyph_count(lit_pd, cl[:, :3], int(a + (b - a) * ease(t))))
        step(int(0.9 * FPS), lambda t: None)
    set_glyph_count(lit_pd, cl[kept, :3], n_keep)
    full_a.VisibilityOff()

    stage.SetInput(f"6  Refit on the kept {keep_pct:.0f} %")
    state_a.VisibilityOn()

    def r6(t):
        k = max(2, int(ease(t) * ns))
        set_traj(k)
        state_src.SetCenter(*sm[k - 1])
        detail.SetInput(f"measurement {k} / {ns}   {nc - n_keep} clusters at the far end dropped")

    step(int(4.5 * FPS), r6)
    state_a.VisibilityOff()
    vtx_a.VisibilityOn()
    full_a.VisibilityOn()
    full_a.GetProperty().SetOpacity(0.35)
    vtxF_a.VisibilityOn()

    stage.SetInput("7  Rescued")
    stage.GetTextProperty().SetColor(0.45, 1.0, 0.55)
    detail.SetInput("blue: full length (fails)    red: FIND (adopted)    Spyral: independent estimate from the track curvature")
    result.GetTextProperty().SetColor(0.95, 0.95, 0.95)
    result.GetTextProperty().SetFontFamilyToCourier()
    result.GetTextProperty().SetJustificationToLeft()
    result.SetDisplayPosition(W - 470, H - 30)
    result.SetInput(f"{'':14s}{'full':>9s}{'FIND':>10s}\n"
                    f"{'KE (MeV)':14s}{full['ke']:9.2f}{fit['ke']:10.2f}\n"
                    f"{'chi2/ndf':14s}{full['chi2ndf']:9.1f}{fit['chi2ndf']:10.3f}\n"
                    f"{'vertex r (mm)':14s}{rF:9.1f}{rK:10.1f}\n"
                    f"{'Spyral KE':14s}{'':9s}{keS:10.2f}")
    step(int(6.0 * FPS), lambda t: None)

if not RESCUE:
    stage.SetInput("4  Kalman filter (GENFIT)")
    state_a.VisibilityOn()


    def s4(t):
        k = max(2, int(ease(t) * ns))
        set_traj(k)
        state_src.SetCenter(*sm[k - 1])
        set_glyph_count(lit_pd, cl[:, :3], int(np.searchsorted(cl_idx, k - 1, side="right")))
        detail.SetInput(f"measurement {k} / {ns}   fitting the helix from the vertex outwards")


    step(int(6.5 * FPS), s4)
    state_a.VisibilityOff()
    vtx_a.VisibilityOn()

    stage.SetInput("5  Result")
    detail.SetInput("fitted trajectory and reaction vertex")
    res = (f"proton   KE  {fit['ke']:.2f} MeV\n"
           f"theta_lab  {fit['theta']:.1f} deg\n"
           f"chi2/ndf  {fit['chi2ndf']:.3f}\n"
           f"vertex z  {vtx[2]:.0f} mm")
    step(int(1.0 * FPS), lambda t: None)
    result.SetInput(res)
    step(int(4.5 * FPS), lambda t: None)

# ---------------------------------------------------------------- write
print(f"{len(frames)} frames")
imageio.mimwrite(out_mp4, frames, fps=FPS, quality=8, macro_block_size=8)
print("wrote", out_mp4)
if out_gif:
    from PIL import Image
    small = [Image.fromarray(f).resize((W * 5 // 8, H * 5 // 8), Image.LANCZOS) for f in frames[::2]]
    pal = small[len(small) // 2].quantize(colors=255, method=Image.Quantize.MEDIANCUT)
    q = [s.quantize(palette=pal, dither=Image.Dither.NONE) for s in small]
    q[0].save(out_gif, save_all=True, append_images=q[1:], duration=int(2000 / FPS), loop=0, optimize=True)
    print("wrote", out_gif)
for k, i in [(f"f{int(100 * x):02d}", int(x * (len(frames) - 1))) for x in (0.1, 0.3, 0.45, 0.55, 0.65, 0.75, 0.85, 0.99)]:
    if i < len(frames):
        imageio.imwrite(out_mp4.rsplit(".", 1)[0] + f"_{k}.png", frames[i])
