#!/usr/bin/env python3
"""Verify the corrected production output across runs 128-139.

Reads AtHit::GetPositionCorr() as written by the reconstruction and checks the three
things the correction is supposed to deliver, per run:

  * the beam reconstructs at the tilt angle measured independently with the magnet off
    (6.47 +- 0.04 deg from run 100),
  * the beam projects back through the detector axis at the entrance window (z ~ 1000 mm),
  * the alpha+alpha elastic opening angle still peaks at 90 deg, i.e. the correction has
    not disturbed the calibration it sits on top of.
"""
import math, pickle, sys, warnings
import numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, "/home/yassid/dec2014_calib")
from batch_scan import load_all, analyse
from cluster_tracks import fit_line, vertex_from_tracks

RUNS = list(range(128, 140))
TB0, K, VZ = 42.3, 0.16 * 2.251 * 10, 22.348   # new TBEntrance=320 origin; mm/tb; mm/us
TILT_REF = 6.473


def peak(v, lo=55.0, hi=125.0):
    w = np.asarray(v); w = w[(w > lo) & (w < hi)]
    if len(w) < 15:
        return float("nan"), 0
    m = np.median(w)
    for _ in range(12):
        s = w[np.abs(w - m) < 25]
        if len(s) < 8:
            break
        m = np.median(s)
    return float(m), int(len(w))


def analyse_run(run):
    evs, groups = load_all(f"/home/yassid/dec2014_calib/prod_{run}.csv")
    cands = []
    for h in groups:
        try:
            r = analyse(h)
        except Exception:
            r = None
        if r:
            cands.append(r)
    bxs, bys, zc, ops, vtx = [], [], [], [], []
    for r in cands:
        T = (r["tb"] - TB0) * 0.16
        P = np.column_stack([r["x"], r["y"], (r["tb"] - TB0) * K])
        try:
            c1, d1 = fit_line(P[r["a1"]], r["q"][r["a1"]])
            c2, d2 = fit_line(P[r["a2"]], r["q"][r["a2"]])
            cb, db = fit_line(P[r["beam"]], r["q"][r["beam"]])
        except Exception:
            continue
        v = vertex_from_tracks([(c1, d1), (c2, d2), (cb, db)])
        if not (0 < v[2] < 1400):
            continue
        vtx.append(v)
        # Regress x,y on drift time rather than a 3D total-least-squares fit: the pads are
        # 8 x 12 mm while the time coordinate is fine, so TLS attributes transverse pad
        # scatter to the direction and biases the angle up (6.89 vs 6.40 deg on run 128).
        mb = r["beam"]
        tbb = r["tb"][mb].astype(float)
        if np.ptp(tbb) > 10 and mb.sum() >= 5:
            bx = np.polyfit((tbb - TB0) * 0.16, r["x"][mb], 1, w=r["q"][mb])[0]
            by = np.polyfit((tbb - TB0) * 0.16, r["y"][mb], 1, w=r["q"][mb])[0]
            bxs.append(bx); bys.append(by)
        # where the beam crosses the detector axis
        m = r["beam"]
        tb = r["tb"][m].astype(float)
        if np.ptp(tb) > 10 and m.sum() >= 5:
            ax = np.polyfit(tb, r["x"][m], 1, w=r["q"][m])
            ay = np.polyfit(tb, r["y"][m], 1, w=r["q"][m])
            den = ax[0] ** 2 + ay[0] ** 2
            if den > 1e-12:
                t = -(ax[0] * ax[1] + ay[0] * ay[1]) / den
                if 0 < t < 600:
                    zc.append((t - TB0) * K)
        # tangent opening angle
        ds = []
        for key, cc in (("a1", c1), ("a2", c2)):
            Q, W = P[r[key]], r["q"][r[key]]
            sel = np.linalg.norm(Q - v, axis=1) < 120.0
            if sel.sum() < 4:
                sel = np.ones(len(Q), bool)
            cx, dx = fit_line(Q[sel], W[sel])
            ds.append(dx if np.dot(dx, cx - v) > 0 else -dx)
        ops.append(math.degrees(math.acos(np.clip(np.dot(ds[0], ds[1]), -1, 1))))
    vtx = np.array(vtx)
    op, nop = peak(ops)
    # Average the slope COMPONENTS, then take the angle once. hypot() is positive definite,
    # so per-event noise inflates every individual angle and a median of angles is biased up
    # (6.96 vs 6.40 deg). The angle of the median slope vector is the unbiased estimate.
    sx, sy = float(np.median(bxs)), float(np.median(bys))
    pol = math.degrees(math.atan2(math.hypot(sx, sy), VZ))
    return dict(run=run, n=len(cands), pol=pol, zc=float(np.median(zc)),
                op=op, nop=nop, vtx=vtx, ops=np.array(ops), sx=sx, sy=sy)


def main():
    res = [analyse_run(r) for r in RUNS]
    print(f"\n{'run':>5} {'cands':>6} {'beam pol':>9} {'z(axis)':>9} {'open pk':>9} {'Nsig':>6}"
          f"   {'vertices in 0<z<1000':>20}")
    for r in res:
        inv = 100 * np.mean((r["vtx"][:, 2] > 0) & (r["vtx"][:, 2] < 1000))
        print(f"{r['run']:>5} {r['n']:>6} {r['pol']:9.2f} {r['zc']:9.0f} {r['op']:9.2f} {r['nop']:6d}"
              f"   {inv:19.0f}%")
    pol = np.array([r["pol"] for r in res]); zc = np.array([r["zc"] for r in res])
    op = np.array([r["op"] for r in res])
    print(f"\n  beam polar   : {pol.mean():.2f} +- {pol.std():.2f} deg   (reference {TILT_REF:.2f} from run 100, B=0)")
    print(f"  axis crossing: {zc.mean():.0f} +- {zc.std():.0f} mm      (entrance window ~1000)")
    print(f"  opening peak : {op.mean():.2f} +- {op.std():.2f} deg   (alpha+alpha elastic = 90)")
    pickle.dump(res, open("/home/yassid/dec2014_calib/prod_results.pkl", "wb"))


if __name__ == "__main__":
    main()
