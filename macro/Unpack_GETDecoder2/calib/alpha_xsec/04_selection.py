#!/usr/bin/env python3
"""Step 4 -- apply the elastic constraint and check that it worked.

The cut is  |theta1 + theta2 - 90| < SUM_WINDOW  (set in config.py).

Makes  out/04_selection.png :
  left   kinematics AFTER the cut -- the band must now sit ON the elastic locus
  right  the two independent determinations of E_cm, which must agree

If the left panel does not tighten onto the white curve, STOP and raise it: something
upstream (drift velocity, entrance, range table) is wrong and nothing downstream is safe.

    python3 04_selection.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config
from lib.loader import build
from lib.eloss import RangeTable

ev = build(config)
tab = RangeTable(config.ELOSS, config.PRESSURE, config.P_TABLE)

keep = np.abs(ev["thsum"] - 90) < config.SUM_WINDOW
if config.MIN_ARM_MM > 0:
    keep &= (ev["L1"] > config.MIN_ARM_MM) & (ev["L2"] > config.MIN_ARM_MM)
print(f"  {keep.sum()} of {len(keep)} events pass ({100*keep.mean():.0f}%)")

th = np.concatenate([ev["th1"][keep], ev["th2"][keep]])
en = np.concatenate([ev["E1"][keep], ev["E2"][keep]])

# E_cm two ways.  They are independent: one uses the arm ranges, the other the depth the
# beam reached before reacting.  Agreement is the check that the energy scale is right.
Ecm_arms = (ev["E1"][keep] + ev["E2"][keep]) / 2.0
Ecm_vtx  = tab.energy_after(ev["path"][keep]) / 2.0
# theta_cm = 2 * theta_lab for equal masses; use the more forward arm
Tcm = 2 * np.minimum(ev["th1"][keep], ev["th2"][keep])

def E_scattered(E, t, m=4.0, M=4.0):
    s = np.sin(np.radians(t)); c = np.cos(np.radians(t))
    disc = M*M - m*m*s*s
    o = np.full(np.shape(t), np.nan); ok = disc >= 0
    o[ok] = E * ((m*c[ok] + np.sqrt(disc[ok]))/(m+M))**2
    return o

grid = np.linspace(0.5, 90, 600)
fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))
ax[0].hist2d(th, np.clip(en, 0, 12), bins=[np.arange(0, 90.5, 0.5), np.arange(0, 12, 0.06)],
             cmap="viridis", cmin=1)
ax[0].plot(grid, E_scattered(config.E_BEAM, grid), "w-", lw=2.2, label="elastic locus")
ax[0].set_xlabel(r"lab angle  $\theta$ [deg]"); ax[0].set_ylabel("energy from range [MeV]")
ax[0].set_title(f"after the cut ({len(th)} tracks)"); ax[0].legend(fontsize=8); ax[0].grid(alpha=.2)

b = np.arange(0, config.E_BEAM/2 + 0.6, 0.1)
ax[1].hist(Ecm_arms, bins=b, color="#4C72B0", label="from arm ranges")
ax[1].hist(Ecm_vtx,  bins=b, color="#DD8452", alpha=.6, label="from vertex depth")
ax[1].axvline(config.E_BEAM/2, color="k", ls="--", lw=1.2,
              label=f"{config.E_BEAM/2:.2f} MeV = $E_b$/2")
ax[1].set_xlabel("$E_{cm}$ [MeV]"); ax[1].set_ylabel("events / 100 keV")
ax[1].set_title("two independent energy scales"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.25)
fig.tight_layout()
out = os.path.join(config.OUT_DIR, "04_selection.png"); fig.savefig(out, dpi=140)

np.savez_compressed(os.path.join(config.OUT_DIR, "selected.npz"),
                    Ecm=Ecm_arms, Ecm_vtx=Ecm_vtx, Tcm=Tcm,
                    th1=ev["th1"][keep], th2=ev["th2"][keep],
                    L1=ev["L1"][keep], L2=ev["L2"][keep], path=ev["path"][keep])
print(f"  E_cm from arm ranges : median {np.median(Ecm_arms):.2f} MeV")
print(f"  E_cm from vertex     : median {np.median(Ecm_vtx):.2f} MeV")
print(f"  above the {config.E_BEAM/2:.2f} MeV ceiling: {100*(Ecm_arms > config.E_BEAM/2).mean():.0f}%")
print(f"  theta_cm: median {np.median(Tcm):.0f} deg, "
      f"{100*((Tcm>70)&(Tcm<110)).mean():.0f}% between 70 and 110")
print(f"  wrote {out} and selected.npz")
