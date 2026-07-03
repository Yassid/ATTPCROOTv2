#!/usr/bin/env python
"""Generate dircluster pseudo-labels on all real events -> real_pseudo.parquet, to test whether
OC can learn real-track separation given MANY (pseudo-)labeled real events (vs 73 gold)."""
import numpy as np, pandas as pd
from dircluster import cluster
R = pd.read_csv("labeling/data/real_events.csv")
rows = []
for ev, g in R.groupby('event'):
    x, y, z, q = (g[c].to_numpy() for c in ['x','y','z','q'])
    if len(x) < 14:
        continue
    lab = cluster(np.stack([x,y,z],1), q, qratio=0.65, min_hits=4)  # geometric pseudo-labels (noise=-1)
    part = np.where(lab < 0, 'noise', 'track')
    for i in range(len(x)):
        rows.append((int(ev), x[i], y[i], z[i], q[i], int(lab[i]), part[i]))
d = pd.DataFrame(rows, columns=['gid','x','y','z','q','label','particle'])
d.to_parquet("data/real_pseudo.parquet", index=False)
print(f"wrote data/real_pseudo.parquet: {d.gid.nunique()} events, {len(d)} hits, "
      f"noise frac {(d.particle=='noise').mean():.2f}")
