#!/usr/bin/env python
"""Sweep smoothed-dE/dx threshold: event-merge% vs per-particle recovery, to find the knee
that keeps merges low AND recovers the 17C recoil. (smooth_q=True.)"""
import sys, numpy as np, pandas as pd
from dircluster import cluster
NEV = int(sys.argv[1]) if len(sys.argv) > 1 else 250
D = pd.read_parquet("data/sim_noisy.parquet")
np.random.seed(0); gids = D.gid.unique(); np.random.shuffle(gids); gids = gids[:NEV]
EV = []
for gid, g in D[D.gid.isin(gids)].groupby('gid'):
    x, y, z, q = (g[c].to_numpy() for c in ['x','y','z','q'])
    part = g['particle'].to_numpy(); tid = g['label'].to_numpy().copy(); tid[part == 'noise'] = -1
    if len(x) < 14: continue
    EV.append((np.stack([x,y,z],1), q, tid, part))

def emerge(pred, tid, mh=8):
    tracks = [t for t in set(tid[tid>=0]) if (tid==t).sum()>=mh]
    for c in set(pred[pred>=0]):
        tt = tid[pred==c]
        if len([t for t in tracks if (tt==t).sum()>=mh]) >= 2: return 1
    return 0

print(f"{'qratio':>6}{'merge%':>8}{'proton%':>9}{'17C%':>7}{'noise_in%':>10}")
for qr in [0.4, 0.5, 0.6, 0.7, 0.8]:
    m = 0; pc=[0,0]; cc=[0,0]; nc=[0,0]
    for P, Q, tid, part in EV:
        pred = cluster(P, Q, qratio=qr, min_hits=4, smooth_q=True)
        m += emerge(pred, tid)
        for pp, acc in [('proton',pc),('17C_recoil',cc),('noise',nc)]:
            msk = part==pp; acc[0]+=int((pred[msk]>=0).sum()); acc[1]+=int(msk.sum())
    print(f"{qr:>6}{100*m/len(EV):>7.1f}%{100*pc[0]/max(pc[1],1):>8.1f}%{100*cc[0]/max(cc[1],1):>6.1f}%"
          f"{100*nc[0]/max(nc[1],1):>9.1f}%")
print("\nwant: low merge%, high proton% & 17C%, low noise_in%")
