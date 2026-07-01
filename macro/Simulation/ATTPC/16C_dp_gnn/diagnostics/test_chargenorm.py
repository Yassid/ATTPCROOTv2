#!/usr/bin/env python
"""Does per-event charge normalization close the sim<->real domain gap?
Baseline raw-charge KS = 0.76. Test schemes that remove absolute gain:
  log      : log10(q)                       (removes multiplicative gain partially)
  zlog     : per-event z-score of log10(q)  (removes gain AND per-event scale)
  rank     : per-event quantile rank in [0,1]
Lower KS = better alignment = more gain-invariant feature for the GNN.
"""
import numpy as np, pandas as pd
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

SC = "/tmp/claude-1000/-home-yassid/5de5b527-5d42-4a6e-95b5-08f19cc1bf97/scratchpad"
real = pd.read_csv(f"{SC}/real_events.csv")
sim  = pd.read_csv(f"{SC}/labeled.csv").rename(columns={"charge": "q"})
# sim: reaction events (proton present), drop unlabeled
sim = sim[sim.trackID >= 0]
pe  = sim.loc[(sim.A == 1) & (sim.Z == 1), "event"].unique()
sim = sim[sim.event.isin(pe)].copy()

# drop nonpositive charge (log safety)
real = real[real.q > 0].copy()
sim  = sim[sim.q  > 0].copy()

def add_norms(df):
    df = df.copy()
    df["log"] = np.log10(df.q)
    g = df.groupby("event")["log"]
    mu = g.transform("mean"); sd = g.transform("std").fillna(0.0)
    df["zlog"] = np.where(sd > 1e-6, (df.log - mu) / sd, 0.0)
    # per-event quantile rank in [0,1]
    df["rank"] = df.groupby("event")["q"].rank(pct=True)
    return df

real = add_norms(real); sim = add_norms(sim)

print(f"{'SCHEME':8s}  {'REAL median [p10,p90]':34s}  {'SIM median [p10,p90]':34s}  KS")
def row(name):
    a, b = real[name].values, sim[name].values
    ks = ks_2samp(a, b).statistic
    print(f"{name:8s}  {np.median(a):7.2f} [{np.percentile(a,10):6.2f},{np.percentile(a,90):7.2f}]"
          f"           {np.median(b):7.2f} [{np.percentile(b,10):6.2f},{np.percentile(b,90):7.2f}]"
          f"           {ks:.3f}")
    return ks
ks = {n: row(n) for n in ["q", "log", "zlog", "rank"]}
print("-"*90)
print(f"baseline raw 'q' KS={ks['q']:.3f}  ->  best normalized KS={min(ks['zlog'],ks['rank']):.3f}")

# plots: overlaid real/sim hist for each scheme
fig, ax = plt.subplots(1, 4, figsize=(19, 4.2))
for i, (name, xlabel) in enumerate([("q","raw charge (log x)"),("log","log10(charge)"),
                                     ("zlog","per-event z(log q)"),("rank","per-event quantile")]):
    a, b = real[name].values, sim[name].values
    if name == "q":
        bins = np.logspace(np.log10(min(a.min(),b.min())), np.log10(max(a.max(),b.max())), 60)
        ax[i].set_xscale("log")
    else:
        bins = np.linspace(min(a.min(),b.min()), max(a.max(),b.max()), 60)
    ax[i].hist(a, bins=bins, density=True, alpha=.5, label="REAL", color="C0")
    ax[i].hist(b, bins=bins, density=True, alpha=.5, label="SIM",  color="C1")
    ax[i].set_title(f"{name}  (KS={ks[name]:.2f})"); ax[i].set_xlabel(xlabel); ax[i].legend()
plt.tight_layout()
out = f"{SC}/chargenorm_test.png"; plt.savefig(out, dpi=110); print(f"saved {out}")
