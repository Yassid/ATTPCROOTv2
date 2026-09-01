#!/usr/bin/env python3
"""Export a handful of events as JSON for the self-contained HTML viewer.

The HTML page (html_viewer_shell.html) has a __DATA__ placeholder; substitute this JSON
into it to get attpc_viewer.html. Kept as two files so the page can be re-generated from a
different run without editing markup:

    python3 export_viewer_json.py hits_run128.csv viewer_events.json
    python3 -c "open('attpc_viewer.html','w').write(
        open('html_viewer_shell.html').read().replace('__DATA__', open('viewer_events.json').read()))"

De-tilting uses the beam axis MEASURED at B=0 (run 100), not the parameter file's
TiltAng/ThetaRot -- see manual section 5.4.8 for why those two differ.
"""
import sys, csv, math, json
from collections import defaultdict
import numpy as np

SRC = sys.argv[1] if len(sys.argv) > 1 else "hits_run128.csv"
DST = sys.argv[2] if len(sys.argv) > 2 else "viewer_events.json"
TILT, AZIM = 5.44, -173.2

t, p = math.radians(TILT), math.radians(AZIM)
b = np.array([math.sin(t) * math.cos(p), math.sin(t) * math.sin(p), math.cos(t)])
zh = np.array([0, 0, 1.0])
v, c = np.cross(b, zh), float(b @ zh)
s = np.linalg.norm(v)
K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
R = np.eye(3) + K + K @ K * ((1 - c) / s**2)

ev = defaultdict(list)
with open(SRC) as fh:
    for r in csv.DictReader(fh):
        ev[int(r["event"])].append((float(r["x"]), float(r["y"]), float(r["z"]),
                                    float(r["xc"]), float(r["yc"]), float(r["q"])))

# A spread of multiplicities: beam-like through reaction-like, so stepping through the
# slider shows real variety rather than ten near-identical beam tracks.
cand = sorted(ev.items(), key=lambda kv: len(kv[1]))
picks, seen = [], set()
for target in (45, 60, 75, 90, 110, 130, 150, 180, 220, 260):
    k, h = min(cand, key=lambda kv: abs(len(kv[1]) - target))
    if k not in seen:
        seen.add(k)
        picks.append((k, h))

out = []
for k, h in picks:
    a = np.array(h)
    cor = np.column_stack([a[:, 3], a[:, 4], a[:, 2]])
    out.append(dict(id=int(k), n=len(h),
                    raw=[[round(x, 1) for x in r] for r in np.column_stack([a[:, 0], a[:, 1], a[:, 2]]).tolist()],
                    cor=[[round(x, 1) for x in r] for r in cor.tolist()],
                    det=[[round(x, 1) for x in r] for r in (cor @ R.T).tolist()],
                    q=[round(x, 1) for x in a[:, 5].tolist()]))

json.dump(dict(run=128, tilt=TILT, azim=AZIM, beam=[round(x, 5) for x in b.tolist()],
               events=out), open(DST, "w"), separators=(",", ":"))
print(f"wrote {DST}: {len(out)} events, {sum(e['n'] for e in out)} hits")
