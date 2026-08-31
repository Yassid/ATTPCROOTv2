import os
#!/usr/bin/env python3
"""Opening-angle distributions vs drift velocity, and vs drift distance.

Left panel: the distribution at several v_D. If the alpha+alpha elastic peak is
real it should be a distinct bump that slides through 90 deg as v_D changes,
with alpha-C / alpha-O elastics sitting at smaller opening angles.

Right panel: at the best v_D, the peak position split by vertex drift distance.
A residual slope there means the Lorentz shear is mis-set, since that is the
one effect that grows with drift time.
"""
import math, pickle, sys, warnings
import numpy as np
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from fit_params import event_geometry, peak_near

cands = pickle.load(open("/home/yassid/dec2014_calib/candidates.pkl", "rb"))
print(f"{len(cands)} candidates")

VS = [1.8, 2.0, 2.2, 2.4, 2.6, 2.8]
fig, ax = plt.subplots(1, 3, figsize=(18, 5.2))

store = {}
for v in VS:
    ops, zs = [], []
    for rec in cands.values():
        try:
            o, vtx, bt = event_geometry(rec, v)
        except Exception:
            continue
        ops.append(o); zs.append(vtx[2])
    store[v] = (np.array(ops), np.array(zs))
    pk, n = peak_near(ops)
    ax[0].hist(ops, bins=np.arange(0, 181, 3), histtype="step", lw=1.6,
               label=f"v_D={v:.1f}  peak={pk:.1f}")

ax[0].axvline(90, color="k", ls="--", lw=1)
ax[0].set_xlabel("opening angle [deg]"); ax[0].set_ylabel("events")
ax[0].set_title("opening angle vs drift velocity")
ax[0].legend(fontsize=8); ax[0].grid(alpha=.3)

# peak position vs v_D
pv = [(v, peak_near(store[v][0])[0]) for v in VS]
ax[1].plot([p[0] for p in pv], [p[1] for p in pv], "o-")
ax[1].axhline(90, color="k", ls="--")
ax[1].set_xlabel("v_D [cm/us]"); ax[1].set_ylabel("peak opening angle [deg]")
ax[1].set_title("peak position vs v_D"); ax[1].grid(alpha=.3)

# peak vs drift distance at the v_D closest to 90
vbest = min(VS, key=lambda v: abs(peak_near(store[v][0])[0] - 90))
ops, zs = store[vbest]
edges = np.percentile(zs, [0, 20, 40, 60, 80, 100])
cx, cy, ce = [], [], []
for lo, hi in zip(edges[:-1], edges[1:]):
    m = (zs >= lo) & (zs < hi)
    pk, n = peak_near(ops[m])
    if not math.isnan(pk):
        cx.append(0.5 * (lo + hi)); cy.append(pk); ce.append(25.0 / max(math.sqrt(n), 1))
ax[2].errorbar(cx, cy, yerr=ce, fmt="o-")
ax[2].axhline(90, color="k", ls="--")
ax[2].set_xlabel("vertex z (drift distance proxy) [mm]")
ax[2].set_ylabel("peak opening angle [deg]")
ax[2].set_title(f"peak vs drift distance at v_D={vbest:.1f}\n(slope => Lorentz shear mis-set)")
ax[2].grid(alpha=.3)

fig.tight_layout()
fig.savefig("/home/yassid/dec2014_calib/opening_dist.png", dpi=110, bbox_inches="tight")
print("wrote opening_dist.png")
print(f"best v_D among scanned: {vbest}")
for v, p in pv:
    print(f"   v_D={v:.2f}  peak={p:.2f}")
