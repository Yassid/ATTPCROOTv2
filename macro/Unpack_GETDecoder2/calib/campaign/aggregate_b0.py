#!/usr/bin/env python3
"""Campaign-wide v_D and tilt from every magnet-OFF run.

At B = 0 the solve is exact (omega*tau = 0, so no E x B drift and no curvature), which is
why the calibration is anchored here rather than on the magnet-on runs. Previously this
came from run 100 alone; now it uses all 68 B = 0 runs, which also tests whether the tilt
really was constant across the campaign -- it was set once mechanically, so any drift with
run number would be a red flag.

The opening angle for equal masses is 90 deg, and the measured peak moves at about
-24 deg per (cm/us), so v_D is where the peak crosses 90. A 5-point grid is enough: the
relation is linear to well under the statistical error over this range.
"""
import glob, os, pickle, sys, math, warnings
import numpy as np
from multiprocessing import Pool
warnings.filterwarnings("ignore")
sys.path.insert(0, "/home/yassid/fair_install/ATTPCROOTv2_fr19port/macro/Unpack_GETDecoder2/calib")
from calibrate_run import geometry, peak

GRID = np.array([2.15, 2.20, 2.25, 2.30, 2.35])
B, E, TILT0 = 0.0, 30000.0, 7.0

def do_run(path):
    run = os.path.basename(path).replace("cands_", "").replace(".pkl", "")
    d = pickle.load(open(path, "rb"))
    cands, nev = d["cands"], d["nevents"]
    ops_g, bts_g = [], []
    for v in GRID:
        ops, bts = [], []
        for rec in cands:
            try:
                o, bt, vtx, _ = geometry(rec, v, B, E, TILT0)
            except Exception:
                continue
            if not (0 < vtx[2] < 1100):
                continue
            ops.append(o); bts.append(bt)
        ops_g.append(np.array(ops)); bts_g.append(np.array(bts))
    return run, nev, len(cands), ops_g, bts_g

def solve(ops_g, bts_g):
    """v_D where the opening-angle peak crosses 90 deg, and the beam tilt there."""
    pks, bts, ns = [], [], []
    for ops, bt in zip(ops_g, bts_g):
        pk, n = peak(ops)
        if math.isnan(pk):
            return None
        sel = (ops > 55) & (ops < 125)
        pks.append(pk); ns.append(n)
        bts.append(float(np.median(bt[sel])) if sel.sum() > 10 else float("nan"))
    A = np.polyfit(GRID, pks, 1)
    v = (90 - A[1]) / A[0]
    return v, float(np.interp(v, GRID, bts)), int(np.mean(ns)), A[0]

if __name__ == "__main__":
    files = sorted(glob.glob("/home/yassid/dec2014_calib/campaign/P300_B0/cands_*.pkl"))
    print(f"aggregating {len(files)} runs", flush=True)
    with Pool(12) as p:
        res = p.map(do_run, files)

    print(f"\n{'run':>5} {'events':>7} {'cands':>6} {'Nsig':>6} {'v_D':>7} {'tilt':>7}")
    per_run, tot_ev, tot_sig = [], 0, 0
    pooled_ops = [[] for _ in GRID]; pooled_bts = [[] for _ in GRID]
    for run, nev, nc, ops_g, bts_g in res:
        tot_ev += nev
        for i in range(len(GRID)):
            pooled_ops[i].append(ops_g[i]); pooled_bts[i].append(bts_g[i])
        s = solve(ops_g, bts_g)
        if s is None:
            print(f"{run:>5} {nev:>7} {nc:>6}   too few"); continue
        v, t, n, _ = s
        tot_sig += n
        per_run.append((run, v, t, n))
        print(f"{run:>5} {nev:>7} {nc:>6} {n:>6} {v:>7.3f} {t:>7.2f}")

    pooled_ops = [np.concatenate(a) for a in pooled_ops]
    pooled_bts = [np.concatenate(a) for a in pooled_bts]
    v, t, n, slope = solve(pooled_ops, pooled_bts)
    vs = np.array([r[1] for r in per_run]); ts = np.array([r[2] for r in per_run])
    print(f"\n=== campaign total: {tot_ev} events, {tot_sig} signal events in {len(per_run)} runs ===")
    print(f"  POOLED   v_D = {v:.4f} cm/us   tilt = {t:.3f} deg   (slope {slope:.1f} deg per cm/us)")
    print(f"  per-run  v_D = {vs.mean():.4f} +- {vs.std(ddof=1)/math.sqrt(len(vs)):.4f} (rms {vs.std(ddof=1):.4f})")
    print(f"  per-run  tilt= {ts.mean():.4f} +- {ts.std(ddof=1)/math.sqrt(len(ts)):.4f} (rms {ts.std(ddof=1):.4f})")
    rn = np.array([int(r[0]) for r in per_run])
    if len(rn) > 3:
        for name, arr in (("v_D", vs), ("tilt", ts)):
            sl, ic = np.polyfit(rn, arr, 1)
            print(f"  drift of {name} with run number: {sl*100:+.4f} per 100 runs")
    pickle.dump(dict(per_run=per_run, pooled=(v, t, n)), open("/home/yassid/dec2014_calib/campaign/b0_result.pkl", "wb"))
