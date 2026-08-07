#!/usr/bin/env python3
"""Fit v_D (and optionally omega*tau, tilt) from the alpha+alpha opening angle.

The 90 deg constraint holds only for elastic scattering off He. Scattering off
the C and O in the CO2 populates smaller opening angles, so we do not average
over everything: we locate the *peak* near 90 deg and ask which v_D puts it
there. Contamination outside the window then shifts nothing.

The degeneracy between v_D, omega*tau and the tilt is broken by the drift
distance: the Lorentz shear grows with drift time while the z-scale error does
not, so requiring the peak to sit at 90 deg in every drift-distance bin
constrains the shear independently of v_D.
"""
import math, pickle, sys, warnings
import numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, "/home/yassid/dec2014_calib")
from cluster_tracks import fit_line, vertex_from_tracks

TB_NS = 160.0
TB0 = 2.2
B_FIELD = 0.5691
TANGENT_R = 120.0     # mm around the vertex used for the tangent direction


def drift_vec(v_d, E, tilt_deg):
    """Langevin components (cm/us). E in V/m; E=None disables the correction."""
    if E is None:
        return 0.0, 0.0, v_d
    t = math.radians(tilt_deg)
    wt = (B_FIELD / E) * (v_d * 1e4)
    f = v_d / (1.0 + wt * wt)
    return (f * wt * math.sin(t),
            f * wt * wt * math.sin(t) * math.cos(t),
            f * (1.0 + wt * wt * math.cos(t) ** 2))


def event_geometry(rec, v_d, E=None, tilt=7.0, tangent_r=TANGENT_R):
    """Opening angle, vertex and beam direction for one candidate event."""
    vx, vy, vz = drift_vec(v_d, E, tilt)
    T = (rec["tb"] - TB0) * TB_NS * 1e-3
    P = np.column_stack([rec["x"] - vx * T * 10.0,
                         rec["y"] - vy * T * 10.0,
                         vz * T * 10.0])
    q = rec["q"]
    segs = []
    for key in ("a1", "a2", "beam"):
        m = rec[key]
        segs.append((m,) + fit_line(P[m], q[m]))
    vtx = vertex_from_tracks([(c, d) for _, c, d in segs])
    dirs = []
    for m, c, d in segs[:2]:
        Q, W = P[m], q[m]
        sel = np.linalg.norm(Q - vtx, axis=1) < tangent_r
        if sel.sum() < 4:
            sel = np.ones(len(Q), dtype=bool)
        c2, d2 = fit_line(Q[sel], W[sel])
        dirs.append(d2 if np.dot(d2, c2 - vtx) > 0 else -d2)
    opening = math.degrees(math.acos(np.clip(np.dot(dirs[0], dirs[1]), -1, 1)))
    db = segs[2][2]
    if db[2] < 0:
        db = -db
    beam_tilt = math.degrees(math.acos(np.clip(abs(db[2]), -1, 1)))
    return opening, vtx, beam_tilt


def peak_near(vals, centre=90.0, half=25.0, lo=55.0, hi=125.0):
    """Robust peak position of the distribution in a window around `centre`."""
    v = np.asarray(vals)
    w = v[(v > lo) & (v < hi)]
    if len(w) < 20:
        return np.nan, 0
    # iterate a trimmed mean toward the mode
    c = np.median(w)
    for _ in range(12):
        sel = w[np.abs(w - c) < half]
        if len(sel) < 10:
            break
        c = np.median(sel)
    return float(c), int(len(w))


def main():
    cands = pickle.load(open("/home/yassid/dec2014_calib/candidates.pkl", "rb"))
    print(f"{len(cands)} candidate events\n")

    print("=== 1D scan: v_D, no Lorentz correction ===")
    print(f"{'v_D':>6} {'peak[deg]':>10} {'N in window':>12} {'median vtx z':>13}")
    rows = []
    for v in np.arange(1.6, 3.21, 0.1):
        ops, zs = [], []
        for rec in cands.values():
            try:
                o, vtx, _ = event_geometry(rec, v)
            except Exception:
                continue
            ops.append(o); zs.append(vtx[2])
        pk, n = peak_near(ops)
        rows.append((v, pk, n))
        print(f"{v:6.2f} {pk:10.2f} {n:12d} {np.median(zs):13.1f}")

    good = [(v, p) for v, p, n in rows if not math.isnan(p)]
    if len(good) >= 2:
        vs = np.array([g[0] for g in good]); ps = np.array([g[1] for g in good])
        i = np.argsort(np.abs(ps - 90))[:2]
        (v1, p1), (v2, p2) = (vs[i[0]], ps[i[0]]), (vs[i[1]], ps[i[1]])
        if p1 != p2:
            vbest = v1 + (90 - p1) * (v2 - v1) / (p2 - p1)
            print(f"\n  -> opening-angle peak reaches 90 deg at v_D = {vbest:.3f} cm/us")


if __name__ == "__main__":
    main()
