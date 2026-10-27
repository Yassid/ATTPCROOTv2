#!/usr/bin/env python3
"""Score ONE track finder against MC truth.  Tracker-agnostic.

Give it a truth file and a candidate pickle -- whatever produced the pickle -- and it runs
the same steps 3-4 the real data go through and reports how often the right answer comes
out.  That is the number 07_corrected_xsec.py divides the cross section by.

  usage: score_tracking.py <truth2.csv> <tracks.pkl> [out.png] [label]

The pickle is {"nevents": N, "cands": [...]}, each candidate a dict with
    tb, x, y, q      equal-length arrays, every hit of the event
    beam, a1, a2     boolean masks into those arrays (beam may be all-False)
    event            the event number, so it can be matched to the truth
which is all a new finder has to produce to be comparable on equal terms.
"""
import sys, os, csv, pickle
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
from lib.geometry import reconstruct
from lib.eloss import RangeTable

truth_csv, pkl = sys.argv[1], sys.argv[2]
out_png = sys.argv[3] if len(sys.argv) > 3 else "sim_efficiency.png"
label = sys.argv[4] if len(sys.argv) > 4 else os.path.basename(pkl)
tab = RangeTable(config.ELOSS, config.PRESSURE, config.P_TABLE)

truth = {}
for r in csv.DictReader(open(truth_csv)):
    truth[int(r["event"])] = (float(r["Ecm"]), float(r["theta_cm"]))

cands = [c for c in pickle.load(open(pkl, "rb"))["cands"] if "event" in c]
if not cands:
    raise SystemExit(f"{pkl} has no candidates carrying an 'event' field -- a finder must "
                     f"keep the event number for truth matching to be possible")

res = {}
for c in cands:
    r = reconstruct(c, config.V_DRIFT, config.T0_BUCKETS, config.TB_NS)
    if r is None:
        continue
    L1, L2 = r["length"]
    res[c["event"]] = dict(
        tcm=2.0 * min(r["theta"]), sum=r["theta_sum"],
        Ecm_arms=(float(tab.energy_of_range(L1)) + float(tab.energy_of_range(L2))) / 2.0,
        Ecm_vtx=float(tab.energy_after(r["path"])) / 2.0,
        hasbeam=bool(c["beam"].sum() >= 2))

evs = sorted(truth)
Et = np.array([truth[e][0] for e in evs]); Tt = np.array([truth[e][1] for e in evs])
has = np.array([e in res for e in evs])
sel = np.array([e in res and abs(res[e]["sum"] - 90.0) < config.SUM_WINDOW for e in evs])
good = np.array([e in res and abs(res[e]["sum"] - 90.0) < config.SUM_WINDOW
                 and abs(res[e]["tcm"] - truth[e][1]) < 10.0 for e in evs])
m = [e for e in evs if e in res]
nan = float("nan")
Er = np.array([res[e]["Ecm_arms"] if e in res else nan for e in evs])
Tr = np.array([res[e]["tcm"] if e in res else nan for e in evs])
dt = np.array([res[e]["tcm"] - truth[e][1] for e in m]) if m else np.array([0.0])
de = np.array([res[e]["Ecm_vtx"] - truth[e][0] for e in m]) if m else np.array([0.0])
da = np.array([res[e]["Ecm_arms"] - truth[e][0] for e in m]) if m else np.array([0.0])
w = lambda v: 0.5 * (np.percentile(v, 84) - np.percentile(v, 16))

print(f"\n  === {label} ===   {len(evs)} simulated reactions")
print(f"    reconstructed             : {100*has.mean():5.1f} %")
print(f"    + 90 deg cut              : {100*sel.mean():5.1f} %")
print(f"    + theta_cm right to 10 deg: {100*good.mean():5.1f} %   <- eps_reco")
if m:
    print(f"    with a beam track         : {100*np.mean([res[e]['hasbeam'] for e in m]):5.1f} %")
    print(f"    theta_cm   residual : {np.median(dt):+6.2f} +- {w(dt):5.2f} deg")
    print(f"    E_cm(arms) residual : {np.median(da):+6.3f} +- {w(da):5.3f} MeV")
    print(f"    E_cm(vtx)  residual : {np.median(de):+6.3f} +- {w(de):5.3f} MeV")

def curve(x, s, edges):
    c = 0.5 * (edges[1:] + edges[:-1])
    num, _ = np.histogram(x[s], bins=edges); den, _ = np.histogram(x, bins=edges)
    eff = np.where(den > 0, num / np.maximum(den, 1), np.nan)
    err = np.sqrt(np.maximum(eff * (1 - eff), 1e-9) / np.maximum(den, 1))
    return c, eff, np.where(den > 0, err, np.nan), den

Eed, Ted = np.arange(0.9, 6.01, 0.3), np.arange(0, 185, 10)
fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.6))
c, eff, err, den = curve(Et, good, Eed)
ax[0].errorbar(c, eff, yerr=err, fmt="o-", ms=3.5, lw=1.3, color="#C44E52")
ax[0].set_xlabel(r"true $E_{cm}$ [MeV]"); ax[0].set_ylabel(r"$\varepsilon_{reco}$")
ax[0].set_ylim(0, 1.05); ax[0].set_title("efficiency vs energy")
c, eff, err, den = curve(Tt, good, Ted)
ax[1].errorbar(c, eff, yerr=err, fmt="o-", ms=3.5, lw=1.3, color="#C44E52")
ax[1].axvspan(80, 100, color="#55A868", alpha=.12, lw=0)
ax[1].set_xlabel(r"true $\theta_{cm}$ [deg]"); ax[1].set_ylabel(r"$\varepsilon_{reco}$")
ax[1].set_ylim(0, 1.05); ax[1].set_title("vs angle (shaded = the band used)")
ax[2].hist(dt, bins=np.arange(-40, 40.5, 2), color="#C44E52", alpha=.85)
ax[2].axvline(0, color="k", lw=1)
ax[2].set_xlabel(r"$\theta_{cm}$ reconstructed $-$ true [deg]"); ax[2].set_ylabel("events")
ax[2].set_title(f"angular resolution ({w(dt):.1f} deg)")
for a in ax: a.grid(alpha=.25)
fig.suptitle(label, fontsize=9)
fig.tight_layout(); fig.savefig(out_png, dpi=140)

np.savez_compressed(os.path.join(os.path.dirname(out_png) or ".", "efficiency.npz"),
                    E=Et, T=Tt, sel=good, has=has, Ereco=Er, Treco=Tr)
band = (Tt > 80) & (Tt < 100)
print(f"\n    in the 80-100 deg band the analysis uses: {int(band.sum())} reactions, "
      f"eps_reco {100*good[band].mean():.1f} %")
print(f"    wrote {out_png} and efficiency.npz")
