#!/usr/bin/env python3
"""Why the trigger-efficiency normalisation is not yet valid.

trigger_eff.py normalises simulation to data on the 1.2-3.0e5 ADC "bulk" and reads the
low-charge deficit as the trigger turn-on. That is only legitimate if the two spectra have
the SAME SHAPE in the normalisation window. They do not, and going from 4k to 40k simulated
events proved the discrepancy is systematic rather than statistical: the ratio in that
window reproduces to two decimals (1.91 / 1.05 / 0.48 at 4k and at 40k) instead of
averaging out.

What the shapes actually are:

  experiment  peaks at 1.8e5, then a SECOND bump near 3.6e5
  simulation  rises smoothly to a hard kinematic edge at ~3.3e5 and stops

The second experimental bump is beam pile-up, and the evidence is quantitative: its charge
is 2.01x the main peak while its hit count is only 1.40x. A second beam particle doubles
the deposited charge but lands on largely the SAME pads, since both follow the beam axis --
so charge doubles and multiplicity does not. The simulation makes one reaction at a time
and cannot produce it: 7.1% of experimental events sit above 4e5 ADC against 0.005% of
simulated ones.

Consequence: the absolute normalisation, and therefore the absolute efficiency scale, is
undetermined. The SHAPE of the low-charge turn-on is still meaningful.

usage: charge_shapes.py [exp.txt] [sim.txt] [out.png]
"""
import sys, os
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

EXP = sys.argv[1] if len(sys.argv) > 1 else "/home/yassid/dec2014_calib/qtot_exp_128.txt"
SIM = sys.argv[2] if len(sys.argv) > 2 else "/home/yassid/dec2014_sim_bulk/qtot_sim_bulk.txt"
OUT = sys.argv[3] if len(sys.argv) > 3 else "/home/yassid/dec2014_calib/plots/charge_shapes.png"


def load(path):
    a = np.loadtxt(path)
    if a.ndim == 1:                       # charge only, no multiplicity column
        return a, None
    return a[:, 0], a[:, 1]


qe, ne = load(EXP)
qs, ns = load(SIM)

fig, ax = plt.subplots(1, 3, figsize=(17, 5))

# --- 1. the two spectra, each normalised to unit area ---------------------------------
a = ax[0]
b = np.logspace(3.6, 6.0, 46)
c = np.sqrt(b[1:] * b[:-1])
he, _ = np.histogram(qe, bins=b)
hs, _ = np.histogram(qs, bins=b)
a.step(c, he / len(qe), where="mid", color="tab:blue", lw=1.8, label=f"experiment ({len(qe)})")
a.step(c, hs / len(qs), where="mid", color="tab:red", lw=1.8, label=f"simulation ({len(qs)})")
a.axvspan(1.2e5, 3.0e5, color="0.85", zorder=0)
a.text(1.85e5, 1.2e-3, "normalisation\nwindow", ha="center", fontsize=8, color="0.35")
a.set_xscale("log"); a.set_yscale("log"); a.set_ylim(1e-4, 0.5)
a.set_xlabel(r"$\Sigma Q_{hit}$ [ADC]"); a.set_ylabel("fraction of events / bin")
a.set_title("The shapes differ inside the window\nused to normalise them")
a.legend(fontsize=8); a.grid(alpha=.3)

# --- 2. linear zoom: the experimental peak and its pile-up replica ---------------------
a = ax[1]
b2 = np.linspace(0, 6e5, 70)
a.hist(qe, bins=b2, histtype="step", color="tab:blue", lw=1.8, density=True, label="experiment")
a.hist(qs, bins=b2, histtype="step", color="tab:red", lw=1.8, density=True, label="simulation")
for x, t in [(1.81e5, "single beam"), (3.62e5, r"$2\times$ = pile-up")]:
    a.axvline(x, color="0.4", ls=":", lw=1.2)
    a.text(x, a.get_ylim()[1] * 0.92, t, rotation=90, fontsize=8, ha="right", va="top", color="0.3")
a.set_xlabel(r"$\Sigma Q_{hit}$ [ADC]"); a.set_ylabel("density")
a.set_title("Experiment has a second bump at exactly $2\\times$\nthe first; simulation has a hard edge")
a.legend(fontsize=8); a.grid(alpha=.3)

# --- 3. the pile-up signature: charge doubles, multiplicity does not -------------------
a = ax[2]
if ne is not None:
    p1 = (qe > 1.60e5) & (qe < 2.00e5)
    p2 = (qe > 3.20e5) & (qe < 4.00e5)
    rq = np.median(qe[p2]) / np.median(qe[p1])
    rn = np.median(ne[p2]) / np.median(ne[p1])
    a.bar([0, 1], [rq, rn], color=["tab:purple", "tab:orange"], width=.55)
    a.axhline(2.0, color="k", ls="--", lw=1.2)
    a.text(1.5, 2.03, "2.0", fontsize=9)
    for i, (v, lab) in enumerate([(rq, f"charge\n{rq:.2f}"), (rn, f"hits\n{rn:.2f}")]):
        a.text(i, v + .05, lab, ha="center", fontsize=10)
    a.set_xticks([0, 1]); a.set_xticklabels(["total charge", "hit multiplicity"])
    a.set_ylim(0, 2.5); a.set_ylabel("second bump / main peak")
    a.set_title("Pile-up signature\ncharge doubles, hits do not: same pads")
    a.grid(alpha=.3, axis="y")

fig.suptitle("Why the efficiency normalisation is not yet valid: the two charge spectra "
             "have different shapes where they are normalised", fontsize=12)
fig.tight_layout()
fig.savefig(OUT, dpi=110, bbox_inches="tight")
print("wrote", OUT)

if ne is not None:
    print(f"pile-up bump / main peak :  charge x{rq:.2f}   hits x{rn:.2f}")
print(f"above 4e5 ADC : experiment {100*(qe>4e5).mean():.2f}%   simulation {100*(qs>4e5).mean():.3f}%")
ok_e, ok_s = ne > 5, ns > 5
if ne is not None and ns is not None:
    print(f"median charge per hit : exp {np.median(qe[ok_e]/ne[ok_e]):.0f}  "
          f"sim {np.median(qs[ok_s]/ns[ok_s]):.0f}  "
          f"ratio {np.median(qs[ok_s]/ns[ok_s])/np.median(qe[ok_e]/ne[ok_e]):.2f}")
