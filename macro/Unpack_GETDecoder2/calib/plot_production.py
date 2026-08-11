#!/usr/bin/env python3
"""Summary figures for the corrected Dec 2014 production (runs 128-139)."""
import math, pickle, sys, warnings
import numpy as np
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

res = pickle.load(open("/home/yassid/dec2014_calib/prod_results.pkl", "rb"))
runs = [r["run"] for r in res]
TILT_REF, TILT_ERR = 6.473, 0.04

fig, ax = plt.subplots(2, 2, figsize=(14, 9))

# 1. beam polar angle per run
a = ax[0][0]
pol = [r["pol"] for r in res]
a.axhspan(TILT_REF - TILT_ERR, TILT_REF + TILT_ERR, color="tab:green", alpha=.25,
          label=f"run 100 (B=0): {TILT_REF:.2f} $\\pm$ {TILT_ERR:.2f}")
a.axhline(TILT_REF, color="tab:green", lw=1.2)
a.plot(runs, pol, "o-", color="tab:blue", label="corrected production")
a.set_xlabel("run"); a.set_ylabel("beam polar angle [deg]")
a.set_title("Beam tilt reproduced run by run\n(reference measured independently with the magnet off)")
a.legend(fontsize=8); a.grid(alpha=.3)

# 2. opening angle per run
a = ax[0][1]
op = [r["op"] for r in res]
a.axhline(90, color="k", ls="--", lw=1.2, label=r"$\alpha+\alpha$ elastic = 90$\degree$")
a.plot(runs, op, "o-", color="tab:red")
a.set_xlabel("run"); a.set_ylabel("opening-angle peak [deg]")
a.set_title("Elastic opening angle after the correction\n(the calibration it sits on top of is undisturbed)")
a.legend(fontsize=8); a.grid(alpha=.3); a.set_ylim(87, 93)

# 3. pooled opening-angle spectrum
a = ax[1][0]
allops = np.concatenate([r["ops"] for r in res])
a.hist(allops, bins=np.arange(0, 181, 3), histtype="stepfilled", color="tab:red", alpha=.35)
a.hist(allops, bins=np.arange(0, 181, 3), histtype="step", color="tab:red", lw=1.4)
a.axvline(90, color="k", ls="--")
a.set_xlabel("opening angle [deg]"); a.set_ylabel("events")
n = sum(r["nop"] for r in res)
a.set_title(f"Pooled opening angle, runs 128-139\n{n} elastic candidates in 55-125$\\degree$; "
            f"peak {np.mean(op):.2f}$\\degree$")
a.grid(alpha=.3)
a.annotate("$\\alpha$+$\\alpha$\nelastic", xy=(90, a.get_ylim()[1]*.55), ha="center", fontsize=9)
a.annotate("pile-up\n(parallel tracks)", xy=(12, a.get_ylim()[1]*.5), ha="center", fontsize=8, color="0.35")
a.annotate("split track", xy=(168, a.get_ylim()[1]*.5), ha="center", fontsize=8, color="0.35")

# 4. vertex distribution
a = ax[1][1]
vtx = np.vstack([r["vtx"] for r in res])
h = a.hist2d(vtx[:, 2], np.hypot(vtx[:, 0], vtx[:, 1]), bins=[np.arange(0, 1400, 25), np.arange(0, 260, 5)],
             cmap="viridis")
a.axvline(1000, color="w", ls="--", lw=1.3)
a.text(1005, 235, "entrance window", color="w", fontsize=8, rotation=90, va="top")
a.set_xlabel("vertex z [mm]"); a.set_ylabel("vertex radius [mm]")
inv = 100 * np.mean((vtx[:, 2] > 0) & (vtx[:, 2] < 1000))
a.set_title(f"Reaction vertices, {len(vtx)} events\n{inv:.0f}% inside the drift volume")
fig.colorbar(h[3], ax=a, label="events")

fig.suptitle("Dec 2014 alphas, runs 128-139 -- reconstruction with the tilt/Lorentz correction applied",
             fontsize=13)
fig.tight_layout()
out = "/home/yassid/dec2014_calib/plots/production_128_139.png"
fig.savefig(out, dpi=110, bbox_inches="tight")
print("wrote", out)
