#!/usr/bin/env python
"""Per-particle hit-level recovery for the splitter: where are track hits lost to gray?"""
import sys, numpy as np, pandas as pd
from dircluster import cluster
NEV = int(sys.argv[1]) if len(sys.argv) > 1 else 250
D = pd.read_parquet("data/sim_noisy.parquet")
np.random.seed(0); gids = D.gid.unique(); np.random.shuffle(gids); gids = gids[:NEV]
agg = {}   # particle -> [clustered_hits, total_hits, ntracks, sizes...]
for gid, g in D[D.gid.isin(gids)].groupby('gid'):
    x, y, z, q = (g[c].to_numpy() for c in ['x','y','z','q'])
    part = g['particle'].to_numpy(); lab_t = g['label'].to_numpy()
    if len(x) < 14: continue
    pred = cluster(np.stack([x,y,z],1), q, qratio=0.4, min_hits=4)
    for p in set(part):
        m = part == p
        a = agg.setdefault(p, [0, 0, [], []])
        a[0] += int((pred[m] >= 0).sum()); a[1] += int(m.sum())
        # per-track sizes for this particle in this event
        for t in set(lab_t[m]):
            tm = m & (lab_t == t); a[2].append(int(tm.sum())); a[3].append(float((pred[tm] >= 0).mean()))
print(f"{'particle':<12}{'hits':>9}{'clustered%':>11}{'ntracks':>9}{'medhits/trk':>12}{'med clustered/trk%':>19}")
for p, (ch, th, sizes, fr) in sorted(agg.items(), key=lambda kv: -kv[1][1]):
    print(f"{p:<12}{th:>9}{100*ch/max(th,1):>10.1f}%{len(sizes):>9}{np.median(sizes) if sizes else 0:>12.0f}"
          f"{100*np.median(fr) if fr else 0:>18.1f}%")
