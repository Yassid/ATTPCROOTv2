#!/usr/bin/env python
"""Package labeled.csv (from noisy-sim digi + extract_labels) into a GNN-ready parquet.

Keeps reaction events (containing the proton ejectile A=1,Z=1). Per-event dense cluster
label; noise hits (trackID=-1) become their OWN label (the "haze" instance, matching the
train_attpc.py convention). Output columns match train_sim.py: gid,x,y,z,q,label,particle,A,Z.
  ~/gnn_env/bin/python package_sim_noisy.py [labeled.csv] [out.parquet]
"""
import sys, numpy as np, pandas as pd

src = sys.argv[1] if len(sys.argv) > 1 else "data/labeled.csv"
out = sys.argv[2] if len(sys.argv) > 2 else "data/sim_noisy.parquet"

d = pd.read_csv(src)
# species per hit: trackID (noise = -1 stays its own species)
d["species"] = d["trackID"].astype(int)

def particle(A, Z):
    if A == 1 and Z == 1: return "proton"
    if Z == 6 and A == 16: return "16C_beam"
    if Z == 6 and A == 17: return "17C_recoil"
    if A < 0: return "noise"
    return f"Z{Z}A{A}"
d["particle"] = [particle(a, z) for a, z in zip(d["A"], d["Z"])]

# keep reaction events (those with a proton hit)
has_p = d.loc[d["particle"] == "proton", "event"].unique()
n_all = d["event"].nunique()
d = d[d["event"].isin(has_p)].copy()
print(f"reaction events (with proton): {d['event'].nunique()}/{n_all}")

# dense per-event cluster label (0..K-1); noise is one of the species -> its own label
d["label"] = d.groupby("event")["species"].transform(lambda s: pd.factorize(s)[0])

d = d.rename(columns={"event": "gid", "charge": "q"})
d = d[["gid", "x", "y", "z", "q", "label", "particle", "A", "Z", "px", "py", "pz"]]
d.to_parquet(out, index=False)

# stats
nev = d["gid"].nunique()
mult = d.groupby("gid")["label"].nunique()
noise_frac = (d["particle"] == "noise").mean()
print(f"wrote {out}: {len(d)} hits, {nev} events, noise hit fraction {noise_frac:.2f}")
print("labels/event distribution:", dict(mult.value_counts().sort_index()))
print("hits by particle:\n", d["particle"].value_counts())
