#!/usr/bin/env python
"""Generic sim-truth quantifier: merge% / recovery / frags-per-track / noise for any
labeling function fn(P,Q)->labels. Compares splitter-only vs splitter->stitcher.
  ~/gnn_env/bin/python quant_eval.py [nEvents]
"""
import sys, numpy as np, pandas as pd
from sklearn.metrics import homogeneity_completeness_v_measure
from dircluster import cluster
from stitch import stitch, attach

NEV = int(sys.argv[1]) if len(sys.argv) > 1 else 300

def metrics(pred, tid, mh=8):
    tracks = [t for t in set(tid[tid >= 0]) if (tid == t).sum() >= mh]
    nmerge = 0
    for c in set(pred[pred >= 0]):
        tt = tid[pred == c]; big = [t for t in tracks if (tt == t).sum() >= mh]
        if len(big) >= 2: nmerge += 1
    frags = []; rec = 0
    for t in tracks:
        tm = tid == t
        cs = [c for c in set(pred[tm & (pred >= 0)]) if ((pred == c) & tm).sum() >= mh]
        frags.append(len(cs))
        clustered = (pred[tm] >= 0).mean()
        merged = any(len([t2 for t2 in tracks if (tid[pred == c] == t2).sum() >= mh]) >= 2 for c in cs)
        if clustered >= 0.6 and not merged: rec += 1
    h, _, _ = homogeneity_completeness_v_measure(tid, pred)
    return nmerge, len(tracks), frags, rec, h

D = pd.read_parquet("data/sim_noisy.parquet")
np.random.seed(0); gids = D.gid.unique(); np.random.shuffle(gids); gids = gids[:NEV]
EV = []
for gid, g in D[D.gid.isin(gids)].groupby('gid'):
    x, y, z, q = (g[c].to_numpy() for c in ['x', 'y', 'z', 'q'])
    tid = g['label'].to_numpy().copy(); tid[g['particle'].to_numpy() == 'noise'] = -1
    if len(x) < 14: continue
    EV.append((np.stack([x, y, z], 1), q, tid))
print(f"events: {len(EV)} (multi-track {sum(len(set(t[t>=0]))>=2 for _,_,t in EV)})\n")

def run(name, fn):
    evm = 0; H = []; F = []; rec = 0; ntr = 0; NC = []
    for P, Q, tid in EV:
        pred = fn(P, Q)
        nmerge, nt, frags, r, h = metrics(pred, tid)
        evm += (nmerge > 0); H.append(h); F += frags; rec += r; ntr += nt
        NC.append(len(set(pred[pred >= 0])))
    print(f"{name:<28} merge {100*evm/len(EV):>5.1f}%  homog {np.mean(H):.3f}  "
          f"frags/trk {np.mean(F):.2f}  recovery {100*rec/max(ntr,1):>5.1f}%  clus/ev {np.mean(NC):.2f}")

print(f"{'method':<28} {'merge':>7} {'homog':>10} {'frags/trk':>10} {'recovery':>10} {'clus/ev':>9}")
run("splitter q=0.4",             lambda P, Q: cluster(P, Q, qratio=0.4, min_hits=4))
run("splitter q=0.4 + ATTACH",    lambda P, Q: attach(P, Q, cluster(P, Q, qratio=0.4, min_hits=4)))
run("splitter q=0.4 + ATTACH+STCH", lambda P, Q: stitch(P, Q, attach(P, Q, cluster(P, Q, qratio=0.4, min_hits=4))))
run("splitter q=0.5 + ATTACH+STCH", lambda P, Q: stitch(P, Q, attach(P, Q, cluster(P, Q, qratio=0.5, min_hits=4))))
print("\ngoal: ATTACH raises recovery (absorbs gray hits); STITCH merges leftover fragments; merge% stays low.")
