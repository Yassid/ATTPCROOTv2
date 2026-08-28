#!/usr/bin/env python3
"""Trigger efficiency from the simulated/experimental total-charge comparison.

The simulation has no trigger, so comparing the total charge per event between simulation
and data isolates it: where the data is deficient relative to simulation, the trigger was
rejecting events. This needs no knowledge of the exclusion-region geometry, which was
never recorded.

Normalisation is taken on the bulk (1.2-3.0e5) rather than on total area, so the
experimental pile-up tail -- events with a second beam particle, which the simulation does
not produce -- cannot distort the low-charge comparison.
"""
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

e = np.loadtxt("qtot_exp_128.txt"); s = np.loadtxt("qtot_sim.txt")
qe, qs = e[:, 0], s[:, 0]
bins = np.logspace(3.8, 6.0, 22); cen = np.sqrt(bins[1:] * bins[:-1])
he, _ = np.histogram(qe, bins=bins); hs, _ = np.histogram(qs, bins=bins)
m = (cen > 1.2e5) & (cen < 3.0e5)
k = (he[m].sum() / len(qe)) / (hs[m].sum() / len(qs))
fe, fs = he / len(qe), hs / len(qs) * k
ee = np.sqrt(np.maximum(he, 1)) / len(qe)
es = np.sqrt(np.maximum(hs, 1)) / len(qs) * k

fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.2))
a = ax[0]
a.step(cen, fe, where="mid", color="tab:blue", lw=1.8, label="experiment (run 128)")
a.step(cen, fs, where="mid", color="tab:red", lw=1.8, label="simulation (no trigger)")
a.set_xscale("log"); a.set_yscale("log")
a.set_xlabel(r"total charge per event  $\Sigma Q_{hit}$  [ADC]")
a.set_ylabel("fraction of events / bin")
a.set_title("Total charge per event\nnormalised on the bulk, not on area")
a.legend(fontsize=8); a.grid(alpha=.3)
a.axvspan(1.2e5, 3.0e5, color="0.85", zorder=0)

a = ax[1]
ok = (fs > 1e-4) & (cen < 4e5)
r = fe[ok] / fs[ok]
re = r * np.sqrt((ee[ok] / np.maximum(fe[ok], 1e-9)) ** 2 + (es[ok] / np.maximum(fs[ok], 1e-9)) ** 2)
a.errorbar(cen[ok], r, yerr=re, fmt="o-", color="tab:green")
a.axhline(1, color="k", ls="--")
a.set_xscale("log"); a.set_ylim(0, 2.2)
a.set_xlabel(r"$\Sigma Q_{hit}$  [ADC]"); a.set_ylabel("experiment / simulation")
a.set_title("Ratio = trigger efficiency\nrises 0.16 -> 1 over $7\\times10^3$ to $6\\times10^4$")
a.grid(alpha=.3)
fig.suptitle("Trigger efficiency from the charge-spectrum comparison, run 128 vs simulation\n"
             "the simulation has no trigger, so the deficit at low charge measures it",
             fontsize=12)
fig.tight_layout()
fig.savefig("/home/yassid/dec2014_calib/plots/trigger_efficiency.png", dpi=110, bbox_inches="tight")
print("wrote plots/trigger_efficiency.png")
print(f"bulk normalisation factor applied to simulation: {k:.3f}")
