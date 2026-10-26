#!/usr/bin/env python3
"""Tracking efficiency against MC truth.

For every simulated reaction we know exactly what was there -- E_cm, theta_cm and both
alpha energies -- so the efficiency is simply: of the events that happened, in how many did
the track finder recover the three segments the analysis needs?

This is the number the cross section has to be divided by, and the one thing the geometric
acceptance cannot supply: eps_geom only asks whether the tracks STAY IN the detector, not
whether the finder SEES them.

  usage: efficiency.py <truth.csv> <tracks.csv> [out.png]
"""
import sys, csv, collections
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

truth_csv, tracks_csv = sys.argv[1], sys.argv[2]
out_png = sys.argv[3] if len(sys.argv) > 3 else "efficiency.png"

truth = {}
for r in csv.DictReader(open(truth_csv)):
    truth[int(r["event"])] = (float(r["Ecm"]), float(r["theta_cm"]),
                              float(r["E1"]), float(r["E2"]))

found = collections.defaultdict(set)
for r in csv.DictReader(open(tracks_csv)):
    found[int(r["event"])].add(int(r["track"]))

ev = sorted(truth)
E   = np.array([truth[e][0] for e in ev])
T   = np.array([truth[e][1] for e in ev])
n   = np.array([len(found.get(e, ())) for e in ev])
ok3 = n >= 3                      # the analysis needs beam + two arms
ok2 = n >= 2
print(f"  {len(ev)} simulated reactions")
print(f"  found >=2 tracks: {100*ok2.mean():.1f}%      >=3 tracks: {100*ok3.mean():.1f}%")

def curve(x, sel, edges):
    c = 0.5 * (edges[1:] + edges[:-1])
    num, _ = np.histogram(x[sel], bins=edges)
    den, _ = np.histogram(x, bins=edges)
    with np.errstate(invalid="ignore", divide="ignore"):
        eff = np.where(den > 0, num / den, np.nan)
        err = np.where(den > 0, np.sqrt(np.maximum(eff * (1 - eff), 1e-9) / np.maximum(den, 1)), np.nan)
    return c, eff, err, den

Eed = np.arange(0.8, 6.0, 0.25)
Ted = np.arange(0, 185, 10)
fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))

for sel, lab, col in ((ok3, r"$\geq$3 tracks (what the analysis needs)", "#C44E52"),
                      (ok2, r"$\geq$2 tracks", "#4C72B0")):
    c, eff, err, den = curve(E, sel, Eed)
    ax[0].errorbar(c, eff, yerr=err, fmt="o-", ms=3.5, lw=1.3, color=col, label=lab)
ax[0].set_xlabel("$E_{cm}$ [MeV]"); ax[0].set_ylabel(r"$\varepsilon_{reco}$")
ax[0].set_ylim(0, 1.05); ax[0].legend(fontsize=8); ax[0].grid(alpha=.25)
ax[0].set_title("tracking efficiency vs energy")

c, eff, err, den = curve(T, ok3, Ted)
ax[1].errorbar(c, eff, yerr=err, fmt="o-", ms=3.5, lw=1.3, color="#C44E52")
ax[1].set_xlabel(r"$\theta_{cm}$ [deg]"); ax[1].set_ylabel(r"$\varepsilon_{reco}$")
ax[1].set_ylim(0, 1.05); ax[1].grid(alpha=.25); ax[1].set_title("vs angle")

h = ax[2].hist2d(E, T, bins=[Eed, Ted], cmap="viridis")
fig.colorbar(h[3], ax=ax[2], label="simulated")
ax[2].set_xlabel("$E_{cm}$ [MeV]"); ax[2].set_ylabel(r"$\theta_{cm}$ [deg]")
ax[2].set_title("where the simulation put events")
fig.tight_layout(); fig.savefig(out_png, dpi=140)
print(f"  wrote {out_png}")

print(f"\n  {'E_cm':>8} {'N':>6} {'eff(>=3)':>9}")
c, eff, err, den = curve(E, ok3, Eed)
for i in range(len(c)):
    if den[i] > 20:
        print(f"  {c[i]:8.2f} {int(den[i]):6d} {eff[i]:9.2f}")
np.savez_compressed(out_png.replace(".png", ".npz"), E=E, T=T, n=n)
