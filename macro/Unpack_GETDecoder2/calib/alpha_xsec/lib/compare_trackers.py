#!/usr/bin/env python3
"""Which track finder is right?  Decided against MC truth, not against each other.

Two finders are available and they disagree about which events are reactions:

  * HDBSCAN in [position, orientation] space      (step 2, pure Python)
  * AtSampleConsensus, ATTPCROOT's RANSAC          (step 2b, needs the framework)

On real data there is no way to tell which is correct -- they simply return different
numbers.  On simulation there is: every event is a known reaction with a known E_cm and
theta_cm, so each finder can be scored on how often the ANALYSIS gets the right answer.

Both are pushed through exactly the same downstream code (geometry.reconstruct and the
90-degree cut), so the comparison is of the finders alone.

  usage: compare_trackers.py <truth2.csv> <ransac_tracks.csv> <hdbscan_tracks.pkl> [out.png]
"""
import sys, os, csv, pickle
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import config
from lib.geometry import reconstruct
from lib.eloss import RangeTable

truth_csv, ransac_csv, hdb_pkl = sys.argv[1], sys.argv[2], sys.argv[3]
out_png = sys.argv[4] if len(sys.argv) > 4 else "compare_trackers.png"
BEAM_ALIGN = 0.85
tab = RangeTable(config.ELOSS, config.PRESSURE, config.P_TABLE)
MM = config.TB_NS * 1e-3 * config.V_DRIFT * 10.0

truth = {}
for r in csv.DictReader(open(truth_csv)):
    truth[int(r["event"])] = (float(r["Ecm"]), float(r["theta_cm"]))

def fit_dir(P, w):
    w = np.asarray(w, float); w = w / w.sum()
    c = (P * w[:, None]).sum(axis=0); D = P - c
    val, vec = np.linalg.eigh((D * w[:, None]).T @ D)
    d = vec[:, np.argmax(val)]
    return c, d / np.linalg.norm(d)

def cands_from_ransac(path):
    ev = {}
    for r in csv.DictReader(open(path)):
        ev.setdefault(int(r["event"]), {}).setdefault(int(r["track"]), []).append(
            (float(r["tb"]), float(r["x"]), float(r["y"]), float(r["q"])))
    out = []
    for e, tracks in ev.items():
        if len(tracks) < 2:
            continue
        segs = []
        for k, rows in tracks.items():
            a = np.asarray(rows)
            P = np.column_stack([a[:, 1], a[:, 2], (a[:, 0] - config.T0_BUCKETS) * MM])
            c, d = fit_dir(P, a[:, 3])
            segs.append(dict(rows=a, d=d, length=float(np.ptp(P @ d)),
                             qmean=float(a[:, 3].mean())))
        beamish = [s for s in segs if abs(s["d"][2]) > BEAM_ALIGN]
        beam = min(beamish, key=lambda s: s["qmean"]) if beamish else None
        arms = sorted([s for s in segs if s is not beam], key=lambda s: -s["length"])[:2]
        if len(arms) < 2:
            continue
        parts, bnd, off, empty = [], [], 0, np.zeros((0, 4))
        for s in (beam, arms[0], arms[1]):
            rows = empty if s is None else s["rows"]
            parts.append(rows); bnd.append((off, off + len(rows))); off += len(rows)
        A = np.vstack(parts)
        def mk(i):
            m = np.zeros(len(A), bool); m[bnd[i][0]:bnd[i][1]] = True; return m
        out.append(dict(event=e, tb=A[:, 0], x=A[:, 1], y=A[:, 2], q=A[:, 3],
                        beam=mk(0), a1=mk(1), a2=mk(2)))
    return out

def cands_from_pickle(path):
    d = pickle.load(open(path, "rb"))
    return [c for c in d["cands"] if "event" in c]

def score(cands, label):
    res = {}
    for c in cands:
        r = reconstruct(c, config.V_DRIFT, config.T0_BUCKETS, config.TB_NS)
        if r is None:
            continue
        L1, L2 = r["length"]
        res[c["event"]] = dict(
            tcm=2.0 * min(r["theta"]), sum=r["theta_sum"], path=r["path"],
            Ecm_arms=(float(tab.energy_of_range(L1)) + float(tab.energy_of_range(L2))) / 2.0,
            Ecm_vtx=float(tab.energy_after(r["path"])) / 2.0,
            hasbeam=bool(c["beam"].sum() >= 2))
    evs = sorted(truth)
    Et = np.array([truth[e][0] for e in evs]); Tt = np.array([truth[e][1] for e in evs])
    has = np.array([e in res for e in evs])
    sel = np.array([e in res and abs(res[e]["sum"] - 90.0) < config.SUM_WINDOW for e in evs])
    m = [e for e in evs if e in res]
    dt = np.array([res[e]["tcm"] - truth[e][1] for e in m]) if m else np.array([0.0])
    de = np.array([res[e]["Ecm_vtx"] - truth[e][0] for e in m]) if m else np.array([0.0])
    da = np.array([res[e]["Ecm_arms"] - truth[e][0] for e in m]) if m else np.array([0.0])
    # how often does the analysis get theta_cm right, not merely return a number?
    good = np.array([e in res and abs(res[e]["sum"] - 90.0) < config.SUM_WINDOW
                     and abs(res[e]["tcm"] - truth[e][1]) < 10.0 for e in evs])
    w = lambda v: 0.5 * (np.percentile(v, 84) - np.percentile(v, 16))
    print(f"\n  === {label} ===")
    print(f"    reconstructed            : {100*has.mean():5.1f} %")
    print(f"    + 90 deg cut             : {100*sel.mean():5.1f} %")
    print(f"    + theta_cm right to 10 d : {100*good.mean():5.1f} %   <- eps_reco")
    print(f"    with a beam track        : {100*np.mean([res[e]['hasbeam'] for e in m]):5.1f} %")
    print(f"    theta_cm  resid : {np.median(dt):+6.2f} +- {w(dt):5.2f} deg")
    print(f"    E_cm(vtx) resid : {np.median(de):+6.3f} +- {w(de):5.3f} MeV")
    print(f"    E_cm(arms)resid : {np.median(da):+6.3f} +- {w(da):5.3f} MeV")
    return dict(E=Et, T=Tt, sel=sel, good=good, has=has, dt=dt, label=label)

A = score(cands_from_ransac(ransac_csv), "AtSampleConsensus (RANSAC)")
B = score(cands_from_pickle(hdb_pkl), "HDBSCAN [position, orientation]")

def curve(x, s, edges):
    c = 0.5 * (edges[1:] + edges[:-1])
    num, _ = np.histogram(x[s], bins=edges); den, _ = np.histogram(x, bins=edges)
    eff = np.where(den > 0, num / np.maximum(den, 1), np.nan)
    err = np.sqrt(np.maximum(eff * (1 - eff), 1e-9) / np.maximum(den, 1))
    return c, eff, np.where(den > 0, err, np.nan)

Eed, Ted = np.arange(0.9, 6.01, 0.3), np.arange(0, 185, 10)
fig, ax = plt.subplots(1, 3, figsize=(16.5, 4.6))
for r, col in ((A, "#C44E52"), (B, "#4C72B0")):
    c, eff, err = curve(r["E"], r["good"], Eed)
    ax[0].errorbar(c, eff, yerr=err, fmt="o-", ms=3.5, lw=1.3, color=col, label=r["label"])
    c, eff, err = curve(r["T"], r["good"], Ted)
    ax[1].errorbar(c, eff, yerr=err, fmt="o-", ms=3.5, lw=1.3, color=col, label=r["label"])
    ax[2].hist(r["dt"], bins=np.arange(-40, 40.5, 2), histtype="step", lw=1.6, color=col,
               label=f"{r['label']}  ({np.median(np.abs(r['dt'])):.1f} deg)")
ax[0].set_xlabel(r"true $E_{cm}$ [MeV]"); ax[0].set_ylabel(r"$\varepsilon_{reco}$")
ax[1].set_xlabel(r"true $\theta_{cm}$ [deg]"); ax[1].set_ylabel(r"$\varepsilon_{reco}$")
ax[2].set_xlabel(r"$\theta_{cm}$ reconstructed $-$ true [deg]"); ax[2].set_ylabel("events")
ax[0].set_title("efficiency vs energy"); ax[1].set_title("efficiency vs angle")
ax[2].set_title("angular resolution"); ax[2].axvline(0, color="k", lw=1)
for a in ax: a.grid(alpha=.25); a.legend(fontsize=7)
ax[0].set_ylim(0, 1.05); ax[1].set_ylim(0, 1.05)
fig.tight_layout(); fig.savefig(out_png, dpi=140)
print(f"\n  wrote {out_png}")
best = A if A["good"].mean() >= B["good"].mean() else B
np.savez_compressed(os.path.join(os.path.dirname(out_png) or ".", "efficiency.npz"),
                    E=best["E"], T=best["T"], sel=best["good"], has=best["has"])
print(f"  efficiency.npz written from: {best['label']}")
