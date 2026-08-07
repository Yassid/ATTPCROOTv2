#!/usr/bin/env python3
"""Full calibration chain for one run: v_D from the 90 deg peak, then the tilt.

For a magnet-OFF run this is exact rather than degenerate:

  * B = 0  =>  omega*tau = 0, so the Langevin drift vector collapses to v_D z-hat.
    There is no E x B lateral drift, hence no dependence on the electric field
    and no coupling to the tilt angle.
  * B = 0  =>  the scattered tracks are straight, so a full-arm line fit is the
    tangent at the vertex. The chord-vs-tangent systematic disappears.

So v_D follows from the alpha+alpha opening angle alone, and the beam direction
then measures the mechanical tilt directly at that v_D.

  usage: calibrate_run.py <hits.root csv dump> [--B <tesla>] [--E <V/m>]
"""
import argparse, csv, math, pickle, sys, warnings
import numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, "/home/yassid/dec2014_calib")
from cluster_tracks import cluster, fit_line, vertex_from_tracks, closest_approach

TB_NS, TB0 = 160.0, 2.2
MIN_ARM_HITS, MIN_ARM_LEN, MAX_VTX_DIST, BEAM_ALIGN = 12, 60.0, 40.0, 0.85


def drift_vec(v_d, B, E, tilt_deg):
    if B == 0.0 or E is None:
        return 0.0, 0.0, v_d
    t = math.radians(tilt_deg)
    wt = (B / E) * (v_d * 1e4)
    f = v_d / (1.0 + wt * wt)
    return (f * wt * math.sin(t), f * wt * wt * math.sin(t) * math.cos(t),
            f * (1.0 + wt * wt * math.cos(t) ** 2))


def load_events(path):
    ev, tb, x, y, q = [], [], [], [], []
    with open(path) as fh:
        for r in csv.reader(fh):
            if r[0] == "event":
                continue
            ev.append(int(r[0])); tb.append(int(r[1]))
            x.append(float(r[2])); y.append(float(r[3])); q.append(float(r[4]))
    ev = np.array(ev); o = np.argsort(ev, kind="stable")
    ev = ev[o]
    data = np.column_stack([np.array(tb)[o], np.array(x)[o], np.array(y)[o], np.array(q)[o]])
    u = np.unique(ev)
    return u, np.split(data, np.searchsorted(ev, u)[1:])


def select(hits, v_ref=2.25):
    tb, x0, y0, q = hits[:, 0], hits[:, 1], hits[:, 2], hits[:, 3]
    P = np.column_stack([x0, y0, (tb - TB0) * TB_NS * 1e-3 * v_ref * 10.0])
    lab, _, _ = cluster(P, q)
    segs = []
    for L in sorted(set(lab) - {-1}):
        m = lab == L
        if m.sum() < 6:
            continue
        c, d = fit_line(P[m], q[m])
        segs.append(dict(mask=m, d=d, c=c, n=int(m.sum()),
                         length=float(np.ptp(P[m] @ d)), qmean=float(q[m].mean())))
    if len(segs) < 3:
        return None
    cand = [s for s in segs if abs(s["d"][2]) > BEAM_ALIGN]
    if not cand:
        return None
    beam = min(cand, key=lambda s: s["qmean"])
    arms = sorted([s for s in segs if s is not beam], key=lambda s: -s["length"])[:2]
    if len(arms) < 2 or any(a["n"] < MIN_ARM_HITS or a["length"] < MIN_ARM_LEN for a in arms):
        return None
    if closest_approach(arms[0]["c"], arms[0]["d"], arms[1]["c"], arms[1]["d"])[1] > MAX_VTX_DIST:
        return None
    return dict(tb=tb, x=x0, y=y0, q=q, beam=beam["mask"], a1=arms[0]["mask"], a2=arms[1]["mask"])


def geometry(rec, v_d, B, E, tilt):
    vx, vy, vz = drift_vec(v_d, B, E, tilt)
    T = (rec["tb"] - TB0) * TB_NS * 1e-3
    P = np.column_stack([rec["x"] - vx * T * 10, rec["y"] - vy * T * 10, vz * T * 10])
    q = rec["q"]
    segs = {k: fit_line(P[rec[k]], q[rec[k]]) for k in ("a1", "a2", "beam")}
    vtx = vertex_from_tracks(list(segs.values()))
    ds = []
    for k in ("a1", "a2"):
        c, d = segs[k]
        ds.append(d if np.dot(d, c - vtx) > 0 else -d)
    op = math.degrees(math.acos(np.clip(np.dot(ds[0], ds[1]), -1, 1)))
    db = segs["beam"][1]
    bt = math.degrees(math.acos(np.clip(abs(db[2]), -1, 1)))
    return op, bt, vtx, db


def peak(vals, lo=55.0, hi=125.0, half=25.0):
    w = np.asarray(vals); w = w[(w > lo) & (w < hi)]
    if len(w) < 20:
        return np.nan, 0
    c = np.median(w)
    for _ in range(12):
        s = w[np.abs(w - c) < half]
        if len(s) < 10:
            break
        c = np.median(s)
    return float(c), int(len(w))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--B", type=float, default=0.0)
    ap.add_argument("--E", type=float, default=30000.0)
    ap.add_argument("--tilt", type=float, default=7.0)
    a = ap.parse_args()

    evs, groups = load_events(a.csv)
    print(f"loaded {len(evs)} events from {a.csv}", flush=True)
    cands = []
    for i, h in enumerate(groups):
        if i % 1000 == 0:
            print(f"  clustering {i}/{len(evs)}  kept={len(cands)}", flush=True)
        try:
            r = select(h)
        except Exception:
            r = None
        if r:
            cands.append(r)
    print(f"selected {len(cands)} candidates\n")

    print(f"B = {a.B} T  ->  Lorentz correction {'DISABLED (exact)' if a.B==0 else 'active'}")
    print(f"{'v_D':>6} {'openpeak':>10} {'Nsig':>6} {'beamtilt':>10}")
    rows = []
    for v in np.arange(1.9, 2.71, 0.05):
        ops, bts = [], []
        for rec in cands:
            try:
                o, bt, vtx, _ = geometry(rec, v, a.B, a.E, a.tilt)
            except Exception:
                continue
            if not (0 < vtx[2] < 1100):
                continue
            ops.append(o); bts.append(bt)
        ops = np.array(ops); bts = np.array(bts)
        pk, n = peak(ops)
        sel = (ops > 55) & (ops < 125)
        bt_med = float(np.median(bts[sel])) if sel.sum() > 10 else float("nan")
        rows.append((v, pk, n, bt_med))
        print(f"{v:6.2f} {pk:10.2f} {n:6d} {bt_med:10.2f}")

    ok = [r for r in rows if not math.isnan(r[1])]
    vs = np.array([r[0] for r in ok]); ps = np.array([r[1] for r in ok]); bs = np.array([r[3] for r in ok])
    A = np.polyfit(vs, ps, 1)
    vbest = (90 - A[1]) / A[0]
    tilt = float(np.interp(vbest, vs, bs))
    print(f"\n  v_D  = {vbest:.3f} cm/us   (peak sensitivity {A[0]:.1f} deg per cm/us)")
    print(f"  tilt = {tilt:.2f} deg   <- beam direction at that v_D"
          f"{'  [exact: B=0, no Lorentz, no curvature]' if a.B==0 else '  [model-dependent]'}")
    with open(a.csv.replace(".csv", "_cands.pkl"), "wb") as fh:
        pickle.dump(cands, fh)


if __name__ == "__main__":
    main()
