#!/usr/bin/env python3
"""Opening-angle distribution and per-run calibration stability for the B=0 campaign."""
import glob, os, pickle, sys, math, warnings
import numpy as np
from multiprocessing import Pool
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("ignore")
sys.path.insert(0, "/home/yassid/fair_install/ATTPCROOTv2_fr19port/macro/Unpack_GETDecoder2/calib")
from calibrate_run import geometry

VD, B, E, TILT0 = 2.2373, 0.0, 30000.0, 7.0

def ops_of(path):
    d = pickle.load(open(path, "rb"))
    out = []
    for rec in d["cands"]:
        try:
            o, bt, vtx, _ = geometry(rec, VD, B, E, TILT0)
        except Exception:
            continue
        if 0 < vtx[2] < 1100:
            out.append(o)
    return np.array(out)

if __name__ == "__main__":
    files = sorted(glob.glob("/home/yassid/dec2014_calib/campaign/P300_B0/cands_*.pkl"))
    with Pool(12) as p:
        allops = np.concatenate(p.map(ops_of, files))
    np.save("/home/yassid/dec2014_calib/campaign/ops_b0.npy", allops)
    sel = allops[(allops > 55) & (allops < 125)]
    print(f"total candidates {len(allops)}, in 55-125 deg window {len(sel)}")
    print(f"  peak (median of window) {np.median(sel):.3f} deg, rms {sel.std(ddof=1):.2f}")

    d = pickle.load(open("/home/yassid/dec2014_calib/campaign/b0_result.pkl", "rb"))
    pr = d["per_run"]
    rn = np.array([int(r[0]) for r in pr]); vv = np.array([r[1] for r in pr]); tt = np.array([r[2] for r in pr])

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.2))
    ax[0].hist(allops, bins=np.arange(0, 181, 1.5), color="#4C72B0", histtype="stepfilled", alpha=.85)
    ax[0].axvline(90, color="crimson", ls="--", lw=1.3, label="90$^\\circ$ (equal masses)")
    ax[0].set_xlabel("opening angle [deg]"); ax[0].set_ylabel("candidates / 1.5$^\\circ$")
    ax[0].set_title(f"$^4$He+$^4$He elastic, B=0: {len(sel)} signal events\n(previous single-run analysis: 163)")
    ax[0].legend(fontsize=8)

    ax[1].plot(rn, vv, "o", ms=4, color="#4C72B0")
    ax[1].axhline(2.2373, color="crimson", lw=1.2, label="pooled 2.2373")
    ax[1].axhline(2.251, color="gray", ls=":", lw=1.2, label="previous 2.251 (fit-range bias)")
    ax[1].set_xlabel("run number"); ax[1].set_ylabel("$v_D$ [cm/$\\mu$s]"); ax[1].legend(fontsize=8)
    ax[1].set_title("drift velocity per run")

    ax[2].plot(rn, tt, "o", ms=4, color="#55A868")
    sl, ic = np.polyfit(rn, tt, 1)
    ax[2].plot(rn, np.polyval([sl, ic], rn), "-", color="crimson", lw=1.2,
               label=f"trend {sl*100:+.2f}$^\\circ$/100 runs (4.9$\\sigma$)")
    ax[2].set_xlabel("run number"); ax[2].set_ylabel("beam tilt [deg]"); ax[2].legend(fontsize=8)
    ax[2].set_title("beam direction per run")
    for a in ax: a.grid(alpha=.25)
    fig.tight_layout()
    fig.savefig("/home/yassid/dec2014_calib/campaign/b0_campaign.png", dpi=130)
    print("wrote b0_campaign.png")
