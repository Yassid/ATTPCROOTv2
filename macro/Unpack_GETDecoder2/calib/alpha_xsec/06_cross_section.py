#!/usr/bin/env python3
"""Step 6 -- the relative differential cross section.

In a thick target the yield in an energy bin is set by how much PATH the beam spends in
that bin, so the stopping power has to be divided back out:

    dsigma/dOmega  proportional to   Y(E, theta) * (dE/dx) / epsilon(E, theta)

This gives a RELATIVE cross section -- the shape in energy, with an arbitrary overall
scale.  Making it absolute needs the integrated beam flux, which is not reliably
available (see the guide).  A correctly acceptance-corrected relative excitation function
is already a result; a wrongly normalised absolute one is not.

Work at theta_cm = 90 deg first.  4He is a spin-0 boson, so for identical particles at
90 deg ALL ODD PARTIAL WAVES VANISH -- the observable is clean and compares directly with
published alpha-alpha phase shifts.  It is also where the acceptance is best.

Makes  out/06_cross_section.png

    python3 06_cross_section.py
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
tab = RangeTable(config.ELOSS, config.PRESSURE, config.P_TABLE)

T_LO, T_HI = 80.0, 100.0           # the theta_cm band around 90 degrees
EDGES = np.arange(1.0, config.E_BEAM / 2 + 0.1, 0.15)
Ec = 0.5 * (EDGES[1:] + EDGES[:-1])

band = (sel["Tcm"] > T_LO) & (sel["Tcm"] < T_HI)
print(f"  {band.sum()} events with theta_cm in {T_LO:.0f}-{T_HI:.0f} deg")

Y, _ = np.histogram(sel["Ecm"][band], bins=EDGES)

# acceptance averaged over the same theta_cm band
jb = (acc["T"] > T_LO) & (acc["T"] < T_HI)
eps_band = acc["eps"][:, jb].mean(axis=1)
eps = np.interp(Ec, acc["E"], eps_band)

dedx = tab.dEdx(2.0 * Ec)          # evaluated at the BEAM energy, 2 * E_cm

good = (eps > 0.05) & (Y > 0)      # do not divide by a tiny acceptance
xs   = np.where(good, Y * dedx / np.where(eps > 0, eps, 1), np.nan)
err  = np.where(good, np.sqrt(np.maximum(Y, 1)) * dedx / np.where(eps > 0, eps, 1), np.nan)

fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
ax[0].step(Ec, Y, where="mid", color="#4C72B0")
ax[0].set_xlabel("$E_{cm}$ [MeV]"); ax[0].set_ylabel("events / 150 keV")
ax[0].set_title(f"raw yield, $\\theta_{{cm}}$ {T_LO:.0f}-{T_HI:.0f}$^\\circ$")

ax[1].plot(Ec, eps, color="#8172B3", lw=1.8)
ax[1].plot(Ec, dedx / dedx.max(), color="#937860", lw=1.4, ls="--",
           label="d$E$/d$x$ (scaled)")
ax[1].set_xlabel("$E_{cm}$ [MeV]"); ax[1].set_ylabel(r"$\varepsilon$")
ax[1].set_ylim(0, 1.05); ax[1].legend(fontsize=8); ax[1].set_title("the two corrections")

ax[2].errorbar(Ec[good], xs[good], yerr=err[good], fmt="o-", ms=3.5, lw=1.2, color="#C44E52")
ax[2].axvline(2.94, color="k", ls="--", lw=1.2, label="$^8$Be 2$^+$ ($E_x$=3.03)")
ax[2].set_xlabel("$E_{cm}$ [MeV]")
ax[2].set_ylabel(r"d$\sigma$/d$\Omega$  [arb.]")
ax[2].set_title(r"relative cross section at $\theta_{cm}\approx90^\circ$")
ax[2].legend(fontsize=8)
for a in ax: a.grid(alpha=.25)
fig.suptitle(f"$E_b$={config.E_BEAM} MeV, {config.PRESSURE} torr, "
             f"$v_D$={config.V_DRIFT} cm/$\\mu$s, 68 runs", fontsize=9, y=1.0)
fig.tight_layout()
out = os.path.join(config.OUT_DIR, "06_cross_section.png"); fig.savefig(out, dpi=140)

print(f"  usable range: E_cm {Ec[good].min():.2f} - {Ec[good].max():.2f} MeV")
print(f"  lowest E_cm reachable is set by the CHAMBER, not by a cut:")
print(f"    beam range {tab.total_mm:.0f} mm vs {config.DRIFT_LENGTH/np.cos(np.radians(config.TILT)):.0f} mm of path available")
print(f"  wrote {out}")
np.savez_compressed(os.path.join(config.OUT_DIR, "cross_section.npz"),
                    E=Ec, yield_=Y, eps=eps, dedx=dedx, xs=xs, err=err)
