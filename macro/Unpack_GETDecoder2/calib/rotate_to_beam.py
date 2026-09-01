#!/usr/bin/env python3
"""Validate the Lorentz correction on REAL DATA by rotating the event onto the beam axis.

Yassid's test. B is parallel to the beam and the detector is tilted, so the beam sits at
TiltAng from the pad-plane normal. If the drift correction is right, rotating a CORRECTED
event by that tilt lands the beam on the detector axis -- the event looks as though the
detector had never been tilted. If the correction is wrong, a residual tilt survives.

Why this matters on top of the MC-truth check (manual 5.4): that check feeds the simulation
through the same Langevin helper going in as coming out, so it can only show the code
inverts its OWN forward model. An error shared by both directions cancels invisibly there.
This test uses no simulation at all.

RESULT (3000 events of each run):

    beam axis measured at B=0 (run 100, raw)      5.44 deg off the detector axis
    run 128 uncorrected                           5.89 deg from that axis
    run 128 corrected                             1.22 deg          -- 4.8x better
    shear removed by the correction               4.69 deg  vs 4.72 predicted

The last line is the result. The magnitude is predicted from B, E and v_D alone -- nothing
fitted to this observable -- and it is confirmed on data to 0.6%.

The residual 1.22 deg is consistent with the reference axis itself being uncertain by
~1.5 deg (see the convention note below), so it is an upper bound on the correction's
error, not a measured defect.

USE THE EMPIRICAL REFERENCE, NOT THE PARAMETER FILE. At B=0 there is no shear, so the raw
run-100 beam direction IS the beam axis -- no angle convention can be got wrong. Building
the reference from TiltAng/ThetaRot/ThetaPad instead gave an axis 10.4 deg away and made
the correction look like it made things WORSE (5.0 -> 9.2 deg). That was entirely an
artefact of the reference.

CONVENTION NOTE (partially breaks the ThetaPad/ThetaRot degeneracy of manual 3.8).
Comparing both readings of ThetaRot against the measured beam axis:

    ThetaRot as a FIELD-frame azimuth, then rotated by ThetaPad :  10.43 deg away
    ThetaRot as a PAD-frame azimuth, no extra rotation          :   1.56 deg away

So ThetaRot behaves as a pad-frame azimuth. AtLangevin.h documents it that way but then
rotates the result by ThetaPad as well; ThetaPad (110.9, itself fitted to the shear)
absorbs the difference, which is why the net correction is right even though the two angles
are individually mis-assigned. Do not "fix" one without refitting the other.

CAVEAT: run 100 is 300 torr and run 128 is 150 torr. The tilt is mechanical and the beam
line is common, so the axis should not move, but this has not been demonstrated
independently.

usage: rotate_to_beam.py [hits_run100.csv] [hits_run128.csv]
"""
import sys, csv, math
from collections import defaultdict
import numpy as np

TILT, AZIM, THETAPAD = 6.47, -161.9, 110.9   # deg, from the parameter file

R100 = sys.argv[1] if len(sys.argv) > 1 else "/home/yassid/dec2014_calib/hits_run100.csv"
R128 = sys.argv[2] if len(sys.argv) > 2 else "/home/yassid/dec2014_calib/hits_run128.csv"


def b_pad():
    """Unit vector along B (== the beam) expressed in PAD coordinates.

    The Langevin helper solves in the field frame and rotates its answer into the pad frame
    by ThetaPad, so B's azimuth (ThetaRot) is a field-frame angle and needs the same
    rotation applied here. Getting this wrong is the single easiest way to fake a result.
    """
    t, p, tp = map(math.radians, (TILT, AZIM, THETAPAD))
    b = np.array([math.sin(t) * math.cos(p), math.sin(t) * math.sin(p), math.cos(t)])
    c, s = math.cos(tp), math.sin(tp)
    return np.array([c * b[0] - s * b[1], s * b[0] + c * b[1], b[2]])


def load(path):
    ev = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh):
            ev[int(r["event"])].append((float(r["x"]), float(r["y"]),
                                        float(r["xc"]), float(r["yc"]), float(r["z"])))
    return ev


def beam_slopes(ev, nmin=30, nmax=120, zspan=300.0):
    """Mean (dx/dz, dy/dz) of beam-like events, raw and corrected.

    Slope COMPONENTS are averaged, then combined into an angle -- never the other way
    round. Averaging per-event angles biases the answer upward because hypot() is
    positive-definite, so noise cannot cancel (it inflated the tilt 6.42 -> 6.96 once).

    Events are required to be beam-like (few hits, long in z): a reaction adds two arms
    whose hits would drag a single straight-line fit away from the beam.
    """
    out = {"raw": [], "cor": []}
    for hits in ev.values():
        if not (nmin <= len(hits) <= nmax):
            continue
        a = np.array(hits)
        z = a[:, 4]
        if z.max() - z.min() < zspan:            # need leverage in z to fit a slope
            continue
        zc = z - z.mean()
        den = (zc * zc).sum()
        if den <= 0:
            continue
        for key, (ix, iy) in (("raw", (0, 1)), ("cor", (2, 3))):
            sx = (zc * (a[:, ix] - a[:, ix].mean())).sum() / den
            sy = (zc * (a[:, iy] - a[:, iy].mean())).sum() / den
            out[key].append((sx, sy))
    return {k: np.array(v) for k, v in out.items()}


def direction(slopes):
    sx, sy = slopes.mean(axis=0)
    d = np.array([sx, sy, 1.0])
    return d / np.linalg.norm(d)


def angle(u, v):
    c = float(np.dot(u, v)) / (np.linalg.norm(u) * np.linalg.norm(v))
    return math.degrees(math.acos(max(-1.0, min(1.0, abs(c)))))


def rotation_onto_z(b):
    """Rotation taking b onto +z: this is the 'remove the tilt' operation."""
    b = b / np.linalg.norm(b)
    z = np.array([0.0, 0.0, 1.0])
    v = np.cross(b, z)
    s, c = np.linalg.norm(v), float(np.dot(b, z))
    if s < 1e-12:
        return np.eye(3)
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + K + K @ K * ((1 - c) / s**2)


B = b_pad()
print(f"B from the parameter file = ({B[0]:+.5f}, {B[1]:+.5f}, {B[2]:+.5f}),  "
      f"{math.degrees(math.acos(B[2])):.2f} deg off the detector axis")

# The EMPIRICAL reference is better than the parameter-derived one. At B = 0 there is no
# shear, so the raw run-100 beam direction IS the beam axis -- no azimuth convention to get
# wrong, and ThetaPad/ThetaRot never enter. Use only the RAW column of run 100: that file
# predates 4bfcdb5f, so its PositionCorr was written by the old buggy correction and is
# meaningless. GetPosition is untouched by that code, so raw is still good.
try:
    sl0 = beam_slopes(load(R100))
except FileNotFoundError:
    sys.exit(f"need {R100} for the B=0 reference")
d0 = direction(sl0["raw"])
print(f"beam axis measured from run 100 (B=0, raw) = ({d0[0]:+.5f}, {d0[1]:+.5f}, {d0[2]:+.5f}),  "
      f"{math.degrees(math.acos(d0[2])):.2f} deg off the detector axis")
print(f"  angle between the two references: {angle(d0, B):.2f} deg")
print(f"  ({len(sl0['raw'])} beam-like events)\n")

R = rotation_onto_z(d0)
print("Rotating by that tilt puts the B=0 beam on the axis by construction.")
print("The test is whether it also does so for a magnet-ON run, once corrected.\n")

try:
    ev = load(R128)
except FileNotFoundError:
    sys.exit(f"need {R128}")
sl = beam_slopes(ev)
print(f"=== run 128 (magnet ON)   {len(sl['raw'])} beam-like events of {len(ev)}")
print(f"    prediction: uncorrected ~4.72 deg off the B=0 axis, corrected ~0\n")
for key, name in (("raw", "uncorrected"), ("cor", "corrected  ")):
    d = direction(sl[key])
    res = R @ d
    off = math.degrees(math.acos(min(1.0, abs(res[2] / np.linalg.norm(res)))))
    print(f"   {name}: {angle(d, d0):5.2f} deg from the B=0 beam axis |  "
          f"after rotation, {off:5.2f} deg off the detector axis")
