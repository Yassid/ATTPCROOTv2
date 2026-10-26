#!/usr/bin/env python3
"""Step 7 -- the cross section with BOTH corrections applied.

Step 6 divides by the geometric acceptance only, which asks whether the two tracks stay
inside the detector.  It does not ask whether the track finder SEES them, and that second
efficiency is the one that has been missing:

    dsigma/dOmega  proportional to  Y * (dE/dx) / (eps_geom * eps_reco)

eps_reco comes from simulation, where the truth is known (lib/efficiency2.py).  Until it is
measured the shape of the excitation function cannot be told apart from the shape of the
reconstruction efficiency -- which matters, because the geometric acceptance turns out to
be FLAT at theta_cm = 90 deg, so geometry is already ruled out as the cause.

    python3 07_corrected_xsec.py            needs out/efficiency.npz from the simulation
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config
from lib.eloss import RangeTable

sel = dict(np.load(os.path.join(config.OUT_DIR, "selected.npz")))
acc = dict(np.load(os.path.join(config.OUT_DIR, "acceptance.npz")))
eff_file = os.path.join(config.OUT_DIR, "efficiency.npz")
if not os.path.exists(eff_file):
    raise SystemExit("no out/efficiency.npz -- run the simulation and lib/efficiency2.py first")
sim = dict(np.load(eff_file))
tab = RangeTable(config.ELOSS, config.PRESSURE, config.P_TABLE)

T_LO, T_HI = 80.0, 100.0
EDGES = np.arange(1.0, config.E_BEAM / 2 + 0.1, 0.2)
Ec = 0.5 * (EDGES[1:] + EDGES[:-1])

# --- measured yield in the 90 deg band ---
band = (sel["Tcm"] > T_LO) & (sel["Tcm"] < T_HI)
Y, _ = np.histogram(sel["Ecm"][band], bins=EDGES)

# --- geometric acceptance over the same band ---
jb = (acc["T"] > T_LO) & (acc["T"] < T_HI)
eps_geo = np.interp(Ec, acc["E"], acc["eps"][:, jb].mean(axis=1))

# --- reconstruction efficiency from the simulation, same band ---
# sim["sel"] is the flag from lib/efficiency2.py: this reaction went through the SAME
# steps 2-4 the data go through and came out passing the 90-degree cut.  Using the pipeline
# result rather than a track count matters -- a count of three tracks can be two arms plus
# a fragment, and says nothing about whether the vertex or the angles came out right.
sb = (sim["T"] > T_LO) & (sim["T"] < T_HI)
num, _ = np.histogram(sim["E"][sb][sim["sel"][sb]], bins=EDGES)
den, _ = np.histogram(sim["E"][sb], bins=EDGES)
with np.errstate(invalid="ignore", divide="ignore"):
    eps_rec = np.where(den > 5, num / den, np.nan)

dedx = tab.dEdx(2.0 * Ec)
tot = eps_geo * eps_rec
good = np.isfinite(tot) & (tot > 0.05) & (Y > 0)
xs = np.where(good, Y * dedx / np.where(good, tot, 1), np.nan)
err = np.where(good, np.sqrt(np.maximum(Y, 1)) * dedx / np.where(good, tot, 1), np.nan)

fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
ax[0].step(Ec, Y, where="mid", color="#4C72B0")
ax[0].set_xlabel("$E_{cm}$ [MeV]"); ax[0].set_ylabel("events"); ax[0].set_title("raw yield")
ax[1].plot(Ec, eps_geo, lw=1.8, color="#8172B3", label=r"$\varepsilon_{geom}$")
ax[1].plot(Ec, eps_rec, lw=1.8, color="#DD8452", label=r"$\varepsilon_{reco}$ (simulation)")
ax[1].plot(Ec, tot, lw=2.2, color="#55A868", label="product")
ax[1].set_xlabel("$E_{cm}$ [MeV]"); ax[1].set_ylim(0, 1.05); ax[1].legend(fontsize=8)
ax[1].set_title("the corrections")
ax[2].errorbar(Ec[good], xs[good], yerr=err[good], fmt="o-", ms=3.5, lw=1.2, color="#C44E52")
ax[2].axvline(2.94, color="k", ls="--", lw=1.2, label="$^8$Be 2$^+$")
ax[2].set_xlabel("$E_{cm}$ [MeV]"); ax[2].set_ylabel(r"d$\sigma$/d$\Omega$ [arb.]")
ax[2].set_title(r"fully corrected, $\theta_{cm}\approx90^\circ$"); ax[2].legend(fontsize=8)
for a in ax: a.grid(alpha=.25)
fig.suptitle(f"$E_b$={config.E_BEAM} MeV, {config.PRESSURE} torr, $v_D$={config.V_DRIFT}", fontsize=9)
fig.tight_layout()
out = os.path.join(config.OUT_DIR, "07_corrected_xsec.png"); fig.savefig(out, dpi=140)
print(f"  usable E_cm {Ec[good].min():.2f} - {Ec[good].max():.2f} MeV")
print(f"  eps_reco spans {np.nanmin(eps_rec):.2f} - {np.nanmax(eps_rec):.2f}")
print(f"  wrote {out}")
