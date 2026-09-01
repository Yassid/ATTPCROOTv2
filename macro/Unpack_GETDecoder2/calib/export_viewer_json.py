#!/usr/bin/env python3
"""Export alpha+alpha candidate events as JSON for the self-contained HTML viewer.

Events are CLUSTERED first and only reaction candidates are kept, because the physics the
viewer is meant to show -- a beam track, a vertex, and two arms leaving it -- is invisible
in an unclustered event picked by hit count alone. The first version of the viewer did
exactly that and showed ten mostly beam-only events with no vertex marked.

For each candidate we export the beam segment, the two arms, the fitted vertex, and the
opening angle between the arms, all in CORRECTED (Lorentz) coordinates and again rotated
onto the detector axis.

What to expect, since it is not what naive intuition suggests:

  * The beam is NOT along z at x=y=0. It is tilted by ~6.4 deg, passing within ~12 mm of
    the detector axis near the entrance and ~95 mm from it by the pad plane.
  * The vertex is NOT at (0,0). Vertices sit a median 50 mm off the axis (90% within
    93 mm), because the beam itself is offset and has a finite spot.
  * The arms are NOT back-to-back. Equal-mass elastic scattering gives a 90 deg opening
    angle in the LAB (measured: 89.3 deg). Back-to-back is the centre-of-mass picture.

usage: export_viewer_json.py [hits_run128.csv] [viewer_events.json] [nEvents]
"""
import sys, csv, math, json
from collections import defaultdict
import numpy as np
sys.path.insert(0, __import__("os").path.dirname(__import__("os").path.abspath(__file__)))
from cluster_tracks import cluster, fit_line, vertex_from_tracks

SRC = sys.argv[1] if len(sys.argv) > 1 else "hits_run128.csv"
DST = sys.argv[2] if len(sys.argv) > 2 else "viewer_events.json"
NWANT = int(sys.argv[3]) if len(sys.argv) > 3 else 10

TILT, AZIM = 5.44, -173.2          # beam axis MEASURED at B=0 (run 100), see manual 5.4.8
t, p = math.radians(TILT), math.radians(AZIM)
bdir = np.array([math.sin(t) * math.cos(p), math.sin(t) * math.sin(p), math.cos(t)])
zh = np.array([0.0, 0.0, 1.0])
v, c = np.cross(bdir, zh), float(bdir @ zh)
s = np.linalg.norm(v)
K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
R = np.eye(3) + K + K @ K * ((1 - c) / s**2)

ev = defaultdict(list)
with open(SRC) as fh:
    for r in csv.DictReader(fh):
        ev[int(r["event"])].append((float(r["xc"]), float(r["yc"]), float(r["z"]), float(r["q"])))

out = []
for k in sorted(ev):
    h = np.array(ev[k])
    if not (60 <= len(h) <= 400):
        continue
    P, q = h[:, :3], h[:, 3]
    try:
        lab, _d, _l = cluster(P, q)
    except Exception:
        continue
    ids = [u for u in set(lab) if u >= 0]
    if len(ids) < 3:
        continue

    segs = []
    for u in ids:
        m = lab == u
        if m.sum() < 8:
            continue
        cc, dd = fit_line(P[m], q[m])
        ext = float(np.ptp(P[m] @ dd))
        segs.append(dict(u=u, m=m, c=cc, d=dd, n=int(m.sum()),
                         qmean=float(q[m].mean()), ext=ext, dz=abs(float(dd[2]))))
    if len(segs) < 3:
        continue

    # beam: best aligned with the drift axis, lowest mean charge among those
    beamc = [t_ for t_ in segs if t_["dz"] > 0.85]
    if not beamc:
        continue
    beam = min(beamc, key=lambda t_: t_["qmean"])
    arms = sorted([t_ for t_ in segs if t_ is not beam], key=lambda t_: -t_["ext"])[:2]
    # BOTH arms must be long enough to see. A 56 mm stub is a valid track but makes the
    # two-arm topology unreadable, which is the whole point of the picture.
    if len(arms) < 2 or arms[1]["ext"] < 120:
        continue

    vtx = vertex_from_tracks([(a["c"], a["d"]) for a in arms] + [(beam["c"], beam["d"])])
    # Vertex inside the active volume: a fit that lands outside it is an extrapolation
    # failure, not a reaction.
    if not (100 < vtx[2] < 1000):
        continue

    # opening angle, with each arm's direction pointed AWAY from the vertex
    ds = []
    for a in arms:
        d = a["d"] if np.dot(a["d"], a["c"] - vtx) > 0 else -a["d"]
        ds.append(d)
    op = math.degrees(math.acos(float(np.clip(np.dot(ds[0], ds[1]), -1, 1))))
    # Keep only the alpha+alpha ELASTIC band. The 90 deg constraint holds for equal-mass
    # scattering off He; scattering off the C and O in the CO2 populates smaller opening
    # angles (60-70 deg is common) and would make the viewer's topology confusing.
    if not (80 < op < 100):
        continue

    def pack(mask):
        return [[round(float(x), 1) for x in row] for row in P[mask].tolist()]

    def packrot(mask):
        return [[round(float(x), 1) for x in row] for row in (P[mask] @ R.T).tolist()]

    out.append(dict(
        id=int(k), n=int(len(h)), op=round(op, 1),
        beam=pack(beam["m"]), a1=pack(arms[0]["m"]), a2=pack(arms[1]["m"]),
        beamR=packrot(beam["m"]), a1R=packrot(arms[0]["m"]), a2R=packrot(arms[1]["m"]),
        vtx=[round(float(x), 1) for x in vtx],
        vtxR=[round(float(x), 1) for x in (R @ vtx)],
        qb=[round(float(x), 1) for x in q[beam["m"]].tolist()],
        q1=[round(float(x), 1) for x in q[arms[0]["m"]].tolist()],
        q2=[round(float(x), 1) for x in q[arms[1]["m"]].tolist()]))
    print(f"  event {k:5d}  {len(h):4d} hits  opening {op:6.2f} deg  "
          f"vertex ({vtx[0]:6.1f},{vtx[1]:6.1f},{vtx[2]:6.1f})  "
          f"arms reach {arms[0]['ext']:5.0f} / {arms[1]['ext']:5.0f} mm")
    if len(out) >= NWANT:
        break

if not out:
    sys.exit("no candidates found")
json.dump(dict(run=128, tilt=TILT, azim=AZIM, beam=[round(float(x), 5) for x in bdir],
               events=out), open(DST, "w"), separators=(",", ":"))
print(f"wrote {DST}: {len(out)} candidates, mean opening "
      f"{np.mean([e['op'] for e in out]):.2f} deg")
