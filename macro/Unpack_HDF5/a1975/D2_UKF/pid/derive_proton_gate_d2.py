#!/usr/bin/env python3
"""Auto-derive a FIRST-PASS proton PID gate for the a1975 D2 16C(d,p)17C channel.

Input : /tmp/pid_d2.txt  (dumped by dumpPID.C: sqrtdEdx brho chi2ndf KE vertexR polar_deg direction)
Output: proton_band_d2.json  (sqrtdEdx-brho polygon, same schema as the proton-target proton_band.json)
        plots/proton_gate_d2.png  (PID plane + clean protons + derived gate)

Method: the proton ejectile of (d,p) is the population the PROTON-hypothesis genfit fits
well -> select "clean protons" (physical KE, on-axis vertex, sane chi2), which are ~99%
backward. Bin them in sqrt(dEdx), take a robust brho envelope (10-90 pct) per populated
bin, pad, and stitch upper+lower edges into a closed polygon. Refine later from full
stats / interactive draw.
"""
import json, os, sys
import numpy as np

SRC = sys.argv[1] if len(sys.argv) > 1 else "/tmp/pid_d2.txt"
HERE = os.path.dirname(os.path.abspath(__file__))
OUT_JSON = os.path.join(HERE, "proton_band_d2.json")
OUT_PNG = os.path.join(HERE, "plots", "proton_gate_d2.png")

d = np.genfromtxt(SRC, names=True)
# clean-proton training sample: physical proton KE, vertex near axis, finite sane chi2/ndf
clean = (d["KE"] > 0.3) & (d["KE"] < 16) & (d["vertexR"] < 60) & (d["chi2ndf"] > 0) & (d["chi2ndf"] < 5)
c = d[clean]
print(f"valid={len(d)}  clean_protons={len(c)}  backward={(c['polar_deg']>=90).sum()}")

x = c["sqrtdEdx"]; y = c["brho"]
# bin in sqrt(dEdx); robust brho envelope per bin
xlo, xhi = np.percentile(x, 1), np.percentile(x, 99)
nb = 16
edges = np.linspace(xlo, xhi, nb + 1)
cx, lo, hi = [], [], []
for i in range(nb):
    m = (x >= edges[i]) & (x < edges[i + 1])
    if m.sum() < 5:           # need enough protons to trust the slice
        continue
    yb = y[m]
    p10, p90 = np.percentile(yb, 10), np.percentile(yb, 90)
    pad = 0.10 * (p90 - p10) + 0.02   # widen a touch so the first-pass gate isn't too tight
    cx.append(0.5 * (edges[i] + edges[i + 1]))
    lo.append(max(0.0, p10 - pad))
    hi.append(p90 + pad)
cx, lo, hi = map(np.array, (cx, lo, hi))

# smooth the envelope edges (3-point moving average) to kill single-bin outlier spikes
def smooth(a):
    if len(a) < 3:
        return a
    k = np.array([0.25, 0.5, 0.25])
    return np.convolve(np.pad(a, 1, mode="edge"), k, mode="valid")
lo, hi = smooth(lo), smooth(hi)
hi = np.maximum(hi, lo + 0.03)  # keep a minimum band width after smoothing
print(f"band bins used: {len(cx)}  sqrtdEdx {cx.min():.1f}-{cx.max():.1f}  brho {lo.min():.3f}-{hi.max():.3f}")

# closed polygon: upper edge low->high sqrtdEdx, then lower edge back
verts = [[float(a), float(b)] for a, b in zip(cx, hi)] + [[float(a), float(b)] for a, b in zip(cx[::-1], lo[::-1])]
gate = {"name": "proton_band_d2", "xaxis": "sqrtdEdx", "yaxis": "brho", "vertices": verts, "Z": 1, "A": 1}
with open(OUT_JSON, "w") as f:
    json.dump(gate, f, indent=4)
print("wrote", OUT_JSON, f"({len(verts)} vertices)")

# validate: fraction of clean protons inside, and of all valid
from matplotlib.path import Path
poly = Path(verts)
inside_clean = poly.contains_points(np.column_stack([c["sqrtdEdx"], c["brho"]])).mean()
allpts = np.column_stack([d["sqrtdEdx"], d["brho"]])
inside_all = poly.contains_points(allpts).sum()
print(f"clean protons inside gate: {inside_clean*100:.0f}%   |  all valid inside: {inside_all}/{len(d)}")

# plot
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(8, 6))
ax.hist2d(d["sqrtdEdx"], d["brho"], bins=[80, 80], range=[[0, 50], [0, 2.0]], cmin=1, cmap="Blues")
ax.scatter(c["sqrtdEdx"], c["brho"], s=6, c="orange", alpha=0.5, label="clean protons (genfit)")
px = [v[0] for v in verts] + [verts[0][0]]; py = [v[1] for v in verts] + [verts[0][1]]
ax.plot(px, py, "r-", lw=2, label="derived proton gate")
ax.set_xlabel("sqrt(dEdx)"); ax.set_ylabel("Brho [T*m]")
ax.set_title("a1975 D2 16C(d,p)17C — first-pass proton PID gate (run_0016, 4k evt)")
ax.legend(); ax.set_xlim(0, 50); ax.set_ylim(0, 2.0)
fig.tight_layout(); fig.savefig(OUT_PNG, dpi=110)
print("wrote", OUT_PNG)
