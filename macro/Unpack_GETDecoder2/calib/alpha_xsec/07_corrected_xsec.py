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

# --- measured yield in the 90 deg band, SIDEBAND SUBTRACTED ---
# The 90-degree cut is strong but not perfect, and the background under it is not a
# constant fraction: measured on this data the purity runs from ~57% at E_cm = 2.25 MeV to
# ~92% at 4.75 MeV.  Subtracting a flat background estimated from the shoulders of the
# theta1+theta2 peak therefore changes the SHAPE, not just the normalisation -- it removes
# ~40% at the bottom of the range and ~8% at the top.  Not subtracting it would flatten
# the excitation function artificially.
#
#   signal   window : |theta1 + theta2 - 90| <  SUM_WINDOW
#   sideband window : SUM_WINDOW <= |...| < 3 * SUM_WINDOW   (twice as wide, so halved)
band = (sel["Tcm"] > T_LO) & (sel["Tcm"] < T_HI)
Y, _ = np.histogram(sel["Ecm"][band], bins=EDGES)
if "Ecm_side" in sel and len(sel["Ecm_side"]):
    B, _ = np.histogram(sel["Ecm_side"], bins=EDGES)
    scale = 1.0 / 2.0                      # sideband is 2x the signal window width
    Yraw = Y.copy()
    Y = np.maximum(Y - scale * B, 0.0)
    pur = np.where(Yraw > 0, Y / np.maximum(Yraw, 1), np.nan)
    print(f"  sideband subtraction: purity {np.nanmin(pur[Yraw>20]):.2f} - "
          f"{np.nanmax(pur[Yraw>20]):.2f} over the bins with >20 raw counts")
else:
    Yraw = Y.copy()
    B = np.zeros_like(Y)
    print("  NOTE: no sideband sample in selected.npz -- yield NOT background subtracted")

# --- geometric acceptance over the same band ---
jb = (acc["T"] > T_LO) & (acc["T"] < T_HI)
eps_geo = np.interp(Ec, acc["E"], acc["eps"][:, jb].mean(axis=1))

# --- reconstruction efficiency from the simulation, same band ---
# Binned on the TRUE E_cm, with the numerator a nested subset of the denominator, so the
# ratio is a well-behaved binomial bounded by 1.
#
# A migration-aware bin-by-bin factor (numerator counted in the RECONSTRUCTED variable)
# would be the more correct object in principle -- the data are binned in the
# reconstructed E_cm, and the resolution is 0.31 MeV against 0.2 MeV bins.  It was tried
# and it is UNUSABLE at this simulation statistics: only ~20 generated reactions land in
# each 0.2 MeV bin inside the 90-degree band, and because its numerator and denominator
# count different event sets the ratio is not bounded by 1 and swings 0.59 - 1.64 from bin
# to bin on pure Poisson noise.  That is noise, not migration structure.
#
# Leaving migration uncorrected is defensible HERE for a specific reason: eps_reco is FLAT
# at 0.85 - 1.00 across E_cm 2.2 - 5.5 MeV, and migration between bins of equal efficiency
# has almost no effect.  It would NOT be defensible near the edges, which is one more
# reason to quote only the flat region.  Enlarging the simulation does NOT fix the
# migration factor -- that was tested -- so if it is ever needed, unfold properly with a
# response matrix rather than a bin-by-bin ratio.
#
# The efficiency is evaluated on COARSE bins and interpolated onto the data's finer grid:
# with ~20 reactions per fine bin the statistical scatter would otherwise dominate the
# correction.  Interpolating is safe precisely because the efficiency is flat.
CO = np.arange(1.0, config.E_BEAM / 2 + 0.5, 0.5)
sb = (sim["T"] > T_LO) & (sim["T"] < T_HI)
den_c, _ = np.histogram(sim["E"][sb], bins=CO)
num_c, _ = np.histogram(sim["E"][sb & sim["sel"]], bins=CO)
Cc = 0.5 * (CO[1:] + CO[:-1])
with np.errstate(invalid="ignore", divide="ignore"):
    eps_c = np.where(den_c >= 10, num_c / np.maximum(den_c, 1), np.nan)
ok_c = np.isfinite(eps_c)
eps_rec = np.interp(Ec, Cc[ok_c], eps_c[ok_c], left=np.nan, right=np.nan)
eps_err = np.interp(Ec, Cc[ok_c],
                    np.sqrt(np.maximum(eps_c[ok_c] * (1 - eps_c[ok_c]), 1e-9)
                            / np.maximum(den_c[ok_c], 1)), left=np.nan, right=np.nan)

dedx = tab.dEdx(2.0 * Ec)
tot = eps_geo * eps_rec
good = np.isfinite(tot) & (tot > 0.05) & (Y > 0)
xs = np.where(good, Y * dedx / np.where(good, tot, 1), np.nan)
# Poisson error on a subtracted yield carries BOTH samples
err = np.where(good, np.sqrt(np.maximum(Yraw, 1) + 0.25 * np.maximum(B, 0))
                     * dedx / np.where(good, tot, 1), np.nan)

fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))
ax[0].step(Ec, Yraw, where="mid", color="#B0B7BF", lw=1.2, label="raw")
ax[0].step(Ec, Y, where="mid", color="#4C72B0", lw=1.6, label="background subtracted")
ax[0].legend(fontsize=8)
ax[0].set_xlabel("$E_{cm}$ [MeV]"); ax[0].set_ylabel("events"); ax[0].set_title("yield")
ax[1].plot(Ec, eps_geo, lw=1.8, color="#8172B3", label=r"$\varepsilon_{geom}$")
ax[1].plot(Ec, eps_rec, lw=1.8, color="#DD8452", label=r"$\varepsilon_{reco}$ (simulation)")
ax[1].fill_between(Ec, eps_rec - eps_err, eps_rec + eps_err, color="#DD8452", alpha=.22, lw=0)
ax[1].axvspan(2.2, 5.5, color="#55A868", alpha=.07, lw=0)
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
print(f"  eps_reco spans {np.nanmin(eps_rec):.2f} - {np.nanmax(eps_rec):.2f}"
      f"  (coarse bins, {int(np.nansum(den_c))} simulated reactions in the band)")
flat = good & (Ec > 2.2) & (Ec < 5.5)
if flat.any():
    print(f"  quote the flat region: E_cm {Ec[flat].min():.2f} - {Ec[flat].max():.2f} MeV, "
          f"eps_reco {np.nanmin(eps_rec[flat]):.2f} - {np.nanmax(eps_rec[flat]):.2f}")
print(f"  wrote {out}")
