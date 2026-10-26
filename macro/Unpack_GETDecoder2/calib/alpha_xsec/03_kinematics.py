#!/usr/bin/env python3
"""Step 3 -- look at the raw kinematics before cutting anything.

Makes  out/03_kinematics.png :
  left   every arm track, lab angle vs energy, with the elastic loci drawn on
  right  the theta1 + theta2 distribution, which is the key plot of this analysis

For equal-mass elastic scattering the two LAB angles measured from the beam satisfy
theta1 + theta2 = 90 deg EXACTLY, at every energy.  Anything away from 90 is not an
elastic pair.  You should see a clear peak at 90 sitting on a broad background.

    python3 03_kinematics.py
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import config
from lib.loader import build

ev = build(config)
print(f"  {len(ev['run'])} reconstructed events, v_D = {config.V_DRIFT} cm/us")

th = np.concatenate([ev["th1"], ev["th2"]])
en = np.concatenate([ev["E1"], ev["E2"]])

def E_scattered(E, theta_deg, m=4.0, M=4.0):
    """Lab energy of the scattered projectile at a given lab angle, elastic."""
    s = np.sin(np.radians(theta_deg)); c = np.cos(np.radians(theta_deg))
    disc = M * M - m * m * s * s
    out = np.full(np.shape(theta_deg), np.nan)
    ok = disc >= 0
    out[ok] = E * ((m * c[ok] + np.sqrt(disc[ok])) / (m + M)) ** 2
    return out

grid = np.linspace(0.5, 90, 600)
fig, ax = plt.subplots(1, 2, figsize=(14, 5.2))

ax[0].hist2d(th, np.clip(en, 0, 12), bins=[np.arange(0, 90.5, 0.5), np.arange(0, 12, 0.06)],
             cmap="viridis", cmin=1)
ax[0].plot(grid, E_scattered(config.E_BEAM, grid), "w-", lw=2.2,
           label=f"elastic, $E_b$={config.E_BEAM}")
for f in (0.7, 0.45, 0.25):
    ax[0].plot(grid, E_scattered(config.E_BEAM * f, grid), "w-", lw=1, alpha=.45)
ax[0].set_xlabel(r"lab angle to the beam  $\theta$ [deg]")
ax[0].set_ylabel("energy from range [MeV]")
ax[0].set_title(f"all arm tracks ({len(th)})")
ax[0].legend(fontsize=8); ax[0].grid(alpha=.2)

ax[1].hist(ev["thsum"], bins=np.arange(0, 200, 1.0), color="#4C72B0")
ax[1].axvline(90, color="#C44E52", lw=1.6, label="90$^\\circ$ = elastic")
for s in (-config.SUM_WINDOW, config.SUM_WINDOW):
    ax[1].axvline(90 + s, color="#C44E52", lw=1, ls="--")
ax[1].set_xlabel(r"$\theta_1 + \theta_2$ [deg]"); ax[1].set_ylabel("events / 1$^\\circ$")
ax[1].set_title("the elastic constraint"); ax[1].legend(fontsize=8); ax[1].grid(alpha=.25)

fig.tight_layout()
out = os.path.join(config.OUT_DIR, "03_kinematics.png")
fig.savefig(out, dpi=140)
keep = np.abs(ev["thsum"] - 90) < config.SUM_WINDOW
print(f"  theta1+theta2 median {np.median(ev['thsum']):.1f} deg")
print(f"  within 90 +- {config.SUM_WINDOW:.0f} deg: {keep.sum()} of {len(keep)} "
      f"({100*keep.mean():.0f}%)")
print(f"  wrote {out}")
