#!/usr/bin/env python3
"""How much does the result depend on the unresolved gas density?

The analysis fails an internal consistency test by 0.77 MeV, and the only surviving lead is
that the EFFECTIVE gas density is ~15% below the logged 299.5 torr.  That is not settled --
a pure density factor does not reproduce the arm-length distribution's shape, and the
pressure was confirmed independently.  So rather than pick a value, this quantifies what
the choice COSTS: the same cross section built at several effective densities.

Nothing here depends on the vertex or the entrance model.  The geometry (arm lengths, beam
path) is read from the cache and is density-independent; only the range table moves.

  usage: density_study.py [out.png]
"""
import sys, os
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
from lib import loader
from lib.eloss import RangeTable

out_png = sys.argv[1] if len(sys.argv) > 1 else "density_study.png"
PRESSURES = (250.0, 275.0, 299.5)
COLS = ("#4C72B0", "#55A868", "#C44E52")

ev = loader.build(config, quiet=True)          # geometry only; energies recomputed below
sim = dict(np.load(os.path.join(config.OUT_DIR, "efficiency.npz")))
EDGES = np.arange(1.0, config.E_BEAM / 2 + 0.1, 0.2)
Ec = 0.5 * (EDGES[1:] + EDGES[:-1])
T_LO, T_HI = 80.0, 100.0

fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.8))
summary = []
for P, col in zip(PRESSURES, COLS):
    tab = RangeTable(config.ELOSS, P, config.P_TABLE)
    E1 = np.array([float(tab.energy_of_range(x)) for x in ev["L1"]])
    E2 = np.array([float(tab.energy_of_range(x)) for x in ev["L2"]])
    Ecm = (E1 + E2) / 2.0
    vtx = tab.energy_after(ev["path"]) / 2.0
    Tcm = 2.0 * np.minimum(ev["th1"], ev["th2"])
    keep = (np.abs(ev["thsum"] - 90) < config.SUM_WINDOW) & (Ecm <= config.E_BEAM / 2)
    band = keep & (Tcm > T_LO) & (Tcm < T_HI)
    Y, _ = np.histogram(Ecm[band], bins=EDGES)
    # the two corrections; eps_geom is 1.0 across the flat region, eps_reco from simulation
    sb = (sim["T"] > T_LO) & (sim["T"] < T_HI)
    CO = np.arange(1.0, config.E_BEAM / 2 + 0.5, 0.5)
    den, _ = np.histogram(sim["E"][sb], bins=CO)
    num, _ = np.histogram(sim["E"][sb & sim["sel"]], bins=CO)
    Cc = 0.5 * (CO[1:] + CO[:-1])
    eps = np.where(den >= 10, num / np.maximum(den, 1), np.nan)
    ok = np.isfinite(eps)
    eps_i = np.interp(Ec, Cc[ok], eps[ok], left=np.nan, right=np.nan)
    dedx = tab.dEdx(2.0 * Ec)
    good = np.isfinite(eps_i) & (eps_i > 0.3) & (Y > 0)
    xs = np.where(good, Y * dedx / np.where(good, eps_i, 1), np.nan)
    xs = xs / np.nanmax(xs)                     # shape only: normalisation is out of scope
    gap = np.median((vtx - Ecm)[keep])
    ax[0].step(Ec, Y, where="mid", color=col, lw=1.4, label=f"{P:.0f} torr")
    ax[1].plot(Ec[good], xs[good], "o-", ms=3.4, lw=1.3, color=col,
               label=f"{P:.0f} torr  (gap {gap:+.2f} MeV)")
    # where the rise sits: the half-maximum crossing of the normalised shape
    half = np.nan
    g = np.flatnonzero(good)
    for i in g[:-1]:
        if np.isfinite(xs[i]) and np.isfinite(xs[i+1]) and xs[i] < 0.5 <= xs[i+1]:
            half = Ec[i] + (0.5 - xs[i]) / (xs[i+1] - xs[i]) * (Ec[i+1] - Ec[i]); break
    summary.append((P, int(keep.sum()), float(np.median(Ecm[keep])), gap, half))
    ax[2].axvline(half, color=col, lw=1.6, label=f"{P:.0f}: {half:.2f} MeV")

ax[0].set_xlabel("$E_{cm}$ [MeV]"); ax[0].set_ylabel("events"); ax[0].set_title("yield")
ax[1].set_xlabel("$E_{cm}$ [MeV]"); ax[1].set_ylabel(r"d$\sigma$/d$\Omega$ [normalised]")
ax[1].set_title(r"corrected shape, $\theta_{cm}\approx90^\circ$")
ax[2].axvspan(2.94 - 0.05, 2.94 + 0.05, color="k", alpha=.5)
ax[2].set_xlim(2.0, 5.0); ax[2].set_ylim(0, 1)
ax[2].set_xlabel("$E_{cm}$ [MeV]"); ax[2].set_yticks([])
ax[2].set_title("half-rise position vs the $^8$Be 2$^+$ (black)")
for a in ax: a.grid(alpha=.25); a.legend(fontsize=8)
fig.suptitle("How much does the unresolved gas density move the answer?", fontsize=10)
fig.tight_layout(); fig.savefig(out_png, dpi=140)

print(f"  {'P [torr]':>9} {'N sel':>7} {'median E_cm':>12} {'vtx-arms':>10} {'half-rise':>10}")
for P, n, med, gap, half in summary:
    print(f"  {P:9.1f} {n:7d} {med:12.2f} {gap:+10.3f} {half:10.2f}")
print(f"\n  the simulation says vtx-arms should be +0.512 MeV")
print(f"  the 8Be 2+ sits at E_cm = 2.94 MeV")
print(f"  wrote {out_png}")
