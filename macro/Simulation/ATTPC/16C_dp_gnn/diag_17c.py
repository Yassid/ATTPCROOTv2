#!/usr/bin/env python
"""Understand the 17C-recoil failure: truth (by particle) vs dircluster pred, in x-y AND z-x
(the recoil is forward/beam-like so its extent is likely along z). Also print charge & geometry.
"""
import sys, numpy as np, pandas as pd
from dircluster import cluster
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

D = pd.read_parquet("data/sim_noisy.parquet")
# events with a sizeable 17C recoil
r17 = D[D.particle == "17C_recoil"].groupby("gid").size()
evs = np.array(r17[r17 >= 20].index.to_numpy(), copy=True)
np.random.seed(2); np.random.shuffle(evs); evs = evs[:5]

pcol = {"proton": "tab:blue", "17C_recoil": "tab:red", "16C_beam": "tab:green", "noise": "lightgray"}
fig, axes = plt.subplots(len(evs), 4, figsize=(15, 3.2*len(evs)))
print(f"{'evt':>6}{'17C hits':>9}{'17C medQ':>9}{'p medQ':>8}{'17C r[mm]':>10}{'17C zspan':>10}{'17C clustered%':>15}")
for r, ev in enumerate(evs):
    g = D[D.gid == ev]; x, y, z, q, part = (g[c].to_numpy() for c in ['x','y','z','q','particle'])
    pred = cluster(np.stack([x,y,z],1), q, qratio=0.65, min_hits=4)
    m17 = part == "17C_recoil"; mp = part == "proton"
    r17v = np.hypot(x[m17], y[m17]); zsp = z[m17].max()-z[m17].min() if m17.any() else 0
    print(f"{ev:>6}{m17.sum():>9}{np.median(q[m17]):>9.0f}{np.median(q[mp]) if mp.any() else 0:>8.0f}"
          f"{np.median(r17v):>10.0f}{zsp:>10.0f}{100*(pred[m17]>=0).mean():>14.1f}%")
    # truth xy, truth zx, pred xy, pred zx
    for a, (xx, yy, tag) in zip(axes[r, :2], [(x, y, "x-y"), (z, x, "z-x")]):
        for p in set(part): a.scatter(xx[part==p], yy[part==p], s=6, c=pcol.get(p, 'k'), label=p)
        a.set_title(f"TRUTH evt{ev} {tag}", fontsize=8); a.tick_params(labelsize=6)
    for a, (xx, yy, tag) in zip(axes[r, 2:], [(x, y, "x-y"), (z, x, "z-x")]):
        a.scatter(xx[pred<0], yy[pred<0], s=4, c='lightgray')
        for i2, cl in enumerate(sorted(set(pred[pred>=0]))):
            mm = pred==cl; a.scatter(xx[mm], yy[mm], s=6, color=cm.tab20.colors[i2%20])
        a.set_title(f"PRED evt{ev} {tag}", fontsize=8); a.tick_params(labelsize=6)
axes[0,0].legend(fontsize=6, markerscale=1.5)
plt.tight_layout(); plt.savefig("diagnostics/diag_17c.png", dpi=95); print("wrote diagnostics/diag_17c.png")
