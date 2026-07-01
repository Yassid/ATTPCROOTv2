#!/usr/bin/env python
"""Sim/overlay-vs-real domain-gap diagnostic for 16C(d,p) GNN training data.
Compares per-hit and per-event feature distributions between:
  - REAL  : run_0016 AtEventH full PSA cloud (beam+reaction+noise mixed)
  - SIM   : regenerated labeled signal clouds (proton + 17C, reaction events)
Outputs a multi-panel PNG + a printed summary table with KS distances.
"""
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import ks_2samp

SC = "/tmp/claude-1000/-home-yassid/5de5b527-5d42-4a6e-95b5-08f19cc1bf97/scratchpad"
real = pd.read_csv(f"{SC}/real_events.csv")
sim  = pd.read_csv(f"{SC}/labeled.csv")

# --- sim: keep only reaction events (those containing a proton: A==1,Z==1) ---
sim = sim[sim.trackID >= 0].copy()
prot_events = sim.loc[(sim.A == 1) & (sim.Z == 1), "event"].unique()
sim_rxn = sim[sim.event.isin(prot_events)].copy()
sim_rxn = sim_rxn.rename(columns={"charge": "q"})

def per_event(df):
    g = df.groupby("event")
    r = pd.DataFrame({
        "nhits":  g.size(),
        "zspan":  g.z.max() - g.z.min(),
        "rmax":   g.apply(lambda d: np.sqrt(d.x**2 + d.y**2).max()),
        "qtot":   g.q.sum(),
    })
    return r

ev_real = per_event(real)
ev_sim  = per_event(sim_rxn)

print("="*72)
print(f"REAL events: {ev_real.shape[0]:5d}   hits: {len(real):8d}  ({len(real)/ev_real.shape[0]:.0f} hits/ev)")
print(f"SIM  events: {ev_sim.shape[0]:5d}   hits: {len(sim_rxn):8d}  ({len(sim_rxn)/max(1,ev_sim.shape[0]):.0f} hits/ev)")
# sim class breakdown
n_p   = ((sim_rxn.A==1)&(sim_rxn.Z==1)).sum()
n_17c = ((sim_rxn.A==17)&(sim_rxn.Z==6)).sum()
print(f"SIM  proton hits: {n_p}  ({100*n_p/len(sim_rxn):.1f}%)   17C hits: {n_17c}  ({100*n_17c/len(sim_rxn):.1f}%)")
prot_per_ev = sim_rxn[(sim_rxn.A==1)&(sim_rxn.Z==1)].groupby("event").size()
print(f"SIM  proton hits/event: median {prot_per_ev.median():.0f}, "
      f"10th {prot_per_ev.quantile(.1):.0f}, 90th {prot_per_ev.quantile(.9):.0f}")
print(f"SIM  shared(ambiguous) hits: {100*sim_rxn.shared.mean():.1f}%")
print("="*72)

# --- per-hit charge KS + per-event feature KS ---
def report(name, a, b):
    ks = ks_2samp(a, b).statistic
    print(f"{name:14s}  REAL med {np.median(a):8.1f} [{np.percentile(a,10):7.1f},{np.percentile(a,90):8.1f}]"
          f"   SIM med {np.median(b):8.1f} [{np.percentile(b,10):7.1f},{np.percentile(b,90):8.1f}]"
          f"   KS={ks:.2f}")
    return ks

print(f"{'FEATURE':14s}  {'REAL median [p10,p90]':38s}   {'SIM median [p10,p90]':38s}  GAP")
report("hit charge",  real.q.values,      sim_rxn.q.values)
report("nhits/event", ev_real.nhits.values, ev_sim.nhits.values)
report("zspan/event", ev_real.zspan.values, ev_sim.zspan.values)
report("rmax/event",  ev_real.rmax.values,  ev_sim.rmax.values)
report("qtot/event",  ev_real.qtot.values,  ev_sim.qtot.values)
print("="*72)
print("KS=0 identical, KS=1 fully disjoint. Large KS => domain gap on that feature.")

# --- plots ---
fig, ax = plt.subplots(2, 3, figsize=(15, 9))
def hist2(a, axis, real_v, sim_v, title, xlabel, logx=False, bins=60):
    if logx:
        real_v = real_v[real_v > 0]; sim_v = sim_v[sim_v > 0]
        bins = np.logspace(np.log10(min(real_v.min(), sim_v.min())),
                           np.log10(max(real_v.max(), sim_v.max())), bins)
        axis.set_xscale("log")
    axis.hist(real_v, bins=bins, density=True, alpha=.5, label="REAL", color="C0")
    axis.hist(sim_v,  bins=bins, density=True, alpha=.5, label="SIM",  color="C1")
    axis.set_title(title); axis.set_xlabel(xlabel); axis.legend()

hist2(0, ax[0,0], real.q.values,        sim_rxn.q.values,      "Per-hit charge",   "charge", logx=True)
hist2(0, ax[0,1], ev_real.nhits.values, ev_sim.nhits.values,   "Hits / event",     "n_hits")
hist2(0, ax[0,2], ev_real.zspan.values, ev_sim.zspan.values,   "z-span / event",   "z-span [mm]")
hist2(0, ax[1,0], ev_real.rmax.values,  ev_sim.rmax.values,    "max xy-radius / event", "r_max [mm]")
hist2(0, ax[1,1], ev_real.qtot.values,  ev_sim.qtot.values,    "total charge / event", "q_tot", logx=True)
# sim proton track length
ax[1,2].hist(prot_per_ev.values, bins=40, color="C3", alpha=.7)
ax[1,2].axvline(prot_per_ev.median(), color="k", ls="--", label=f"median {prot_per_ev.median():.0f}")
ax[1,2].set_title("SIM proton hits / event (target class size)")
ax[1,2].set_xlabel("proton n_hits"); ax[1,2].legend()

plt.tight_layout()
out = f"{SC}/dist_compare.png"
plt.savefig(out, dpi=110)
print(f"\nsaved {out}")
