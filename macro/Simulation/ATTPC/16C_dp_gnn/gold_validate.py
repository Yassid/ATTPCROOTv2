#!/usr/bin/env python
"""Validate dircluster (and dircluster->stitch) on the REAL hand-labeled gold set.
Reports merge/recovery/ARI vs manual labels for splitter-only vs splitter+stitcher,
and a GOLD | split | stitched gallery.
"""
import numpy as np, pandas as pd
from sklearn.metrics import adjusted_rand_score, homogeneity_completeness_v_measure
from dircluster import cluster
from stitch import stitch
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import cm

G = pd.read_parquet("labeling/data/labels.parquet"); G = G[G.reviewed]
evs = sorted(G.event.unique()); MH = 8

def metrics(pred, gold):
    tracks = [t for t in set(gold) if (gold == t).sum() >= MH]
    merged = 0
    for c in set(pred[pred >= 0]):
        tt = gold[pred == c]
        if len([t for t in tracks if (tt == t).sum() >= MH]) >= 2: merged = 1
    rec = 0
    for t in tracks:
        tm = gold == t; ok = 0
        for c in set(pred[tm & (pred >= 0)]):
            inter = ((pred == c) & tm).sum()
            if inter / tm.sum() >= 0.6 and inter / (pred == c).sum() >= 0.6: ok = 1
        rec += ok
    return merged, len(tracks), rec

def evalfn(name, fn):
    M = []; ntr = 0; rec = 0; ARI = []; HOM = []
    for ev in evs:
        g = G[G.event == ev]; x, y, z, q = (g[c].to_numpy() for c in ['x','y','z','q'])
        gold = g['label'].to_numpy()
        if len(x) < 14: continue
        pred = fn(np.stack([x, y, z], 1), q)
        m, nt, r = metrics(pred, gold); M.append(m); ntr += nt; rec += r
        ARI.append(adjusted_rand_score(gold, pred)); HOM.append(homogeneity_completeness_v_measure(gold, pred)[0])
    print(f"{name:<24} merge {100*np.mean(M):>4.1f}%   recovery {100*rec/max(ntr,1):>5.1f}% ({rec}/{ntr})   "
          f"ARI {np.mean(ARI):.3f}   homog {np.mean(HOM):.3f}")

split = lambda P, Q: cluster(P, Q, qratio=0.65, min_hits=4)
stch  = lambda P, Q: stitch(P, Q, cluster(P, Q, qratio=0.65, min_hits=4))
print(f"gold reviewed events: {len(evs)}")
evalfn("splitter only", split)
evalfn("splitter + STITCH", stch)

# gallery: GOLD | split | stitched
sel = evs[:8]
fig, axes = plt.subplots(len(sel), 3, figsize=(10.5, 3.0*len(sel)))
for r2, ev in enumerate(sel):
    g = G[G.event == ev]; x, y, z, q = (g[c].to_numpy() for c in ['x','y','z','q'])
    ps = cluster(np.stack([x, y, z], 1), q, qratio=0.65, min_hits=4)
    pt = stitch(np.stack([x, y, z], 1), q, ps)
    for a, lab, tag in [(axes[r2,0], g.label.to_numpy(), "GOLD"), (axes[r2,1], ps, "split"), (axes[r2,2], pt, "stitched")]:
        a.scatter(x[lab < 0], y[lab < 0], s=4, c='lightgray')
        for i, cl in enumerate(sorted(set(lab[lab >= 0]))):
            m = lab == cl; a.scatter(x[m], y[m], s=8, color=cm.tab10.colors[i % 10])
        a.set_title(f"{tag} evt{ev} ({len(set(lab[lab>=0]))} tr)", fontsize=8); a.tick_params(labelsize=6)
plt.tight_layout(); plt.savefig("diagnostics/gold_stitch.png", dpi=95)
print("wrote diagnostics/gold_stitch.png")
