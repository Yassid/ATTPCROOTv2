#!/usr/bin/env python3
"""Drift-velocity calibration from the alpha+alpha elastic opening angle.

For equal-mass non-relativistic elastic scattering the two outgoing particles
always separate by 90 deg. The pad coordinates (x0,y0) are fixed by geometry,
but z is proportional to the drift velocity, so v_D is the one knob that changes
the reconstructed opening angle. We scan it for the value that gives 90 deg.

Segmentation here is deliberately explicit rather than a general clusterer: this
is one hand-picked clean event, and the point is to see the numbers.
"""
import csv, math
import numpy as np

from view_event import drift_vector, TB_NS, TB_ENTRANCE, Z_PADPLANE, TILT, B_FIELD, E_FIELD


def load(event, path="/home/yassid/dec2014_calib/events_run0128.csv"):
    r = [row for row in csv.DictReader(open(path)) if int(row["event"]) == event]
    return (np.array([int(x["tb"]) for x in r]),
            np.array([float(x["x"]) for x in r]),
            np.array([float(x["y"]) for x in r]),
            np.array([float(x["q"]) for x in r]))


def positions(tb, x0, y0, v_d, lorentz=True):
    vx, vy, vz, _ = drift_vector(v_d)
    T = (TB_ENTRANCE - tb.astype(float)) * TB_NS * 1e-3      # us
    if lorentz:
        return (x0 - vx * T * 10.0, y0 - vy * T * 10.0, Z_PADPLANE - vz * T * 10.0)
    return (x0.copy(), y0.copy(), Z_PADPLANE - v_d * T * 10.0)


def fit_line(P, w):
    """Charge-weighted total-least-squares 3D line. Returns (centroid, unit dir)."""
    w = w / w.sum()
    c = (P * w[:, None]).sum(axis=0)
    D = P - c
    # weighted covariance -> principal axis
    C = (D * w[:, None]).T @ D
    val, vec = np.linalg.eigh(C)
    d = vec[:, np.argmax(val)]
    return c, d / np.linalg.norm(d)


def segment(tb, x0, y0, q):
    """Split event 258 into beam / track A / track B by time bucket and pad x."""
    good = tb > 150                       # drop the 3 stray low-tb noise hits
    beam = good & (tb >= 238)
    vtx_tb = 229.0
    outg = good & (tb < vtx_tb)
    # the two arms separate cleanly in x about the vertex x
    A = outg & (x0 < 93.0)
    B = outg & (x0 >= 93.0)
    return beam, A, B


def opening(event=258, lorentz=True, verbose=False):
    tb, x0, y0, q = load(event)
    beam, A, B = segment(tb, x0, y0, q)
    if verbose:
        print(f"  segmentation: beam={beam.sum()}  A={A.sum()}  B={B.sum()}  dropped={(tb<=150).sum()}")

    def ang(v_d):
        X, Y, Z = positions(tb, x0, y0, v_d, lorentz)
        P = np.column_stack([X, Y, Z])
        cA, dA = fit_line(P[A], q[A])
        cB, dB = fit_line(P[B], q[B])
        cb, db = fit_line(P[beam], q[beam])
        # orient both arms away from the vertex (approx: mean of the two centroids)
        vtx = 0.5 * (cA + cB)
        if np.dot(dA, cA - vtx) < 0: dA = -dA
        if np.dot(dB, cB - vtx) < 0: dB = -dB
        if np.dot(db, cb - vtx) > 0: db = -db      # beam points toward the vertex
        return (math.degrees(math.acos(np.clip(np.dot(dA, dB), -1, 1))), dA, dB, db, cA, cB)

    return ang


def main():
    print("=== event 258: opening angle vs drift velocity ===\n")
    for lorentz in (False, True):
        ang = opening(258, lorentz=lorentz, verbose=not lorentz)
        tag = "WITH Lorentz" if lorentz else "no Lorentz  "
        print(f"\n  {tag}")
        print(f"   {'v_D [cm/us]':>12} {'opening [deg]':>14}")
        best, bestd = None, 1e9
        for v in np.arange(1.6, 4.01, 0.05):
            a = ang(v)[0]
            if abs(a - 90) < bestd:
                bestd, best = abs(a - 90), v
            if abs(v * 20 - round(v * 20)) < 1e-9 and abs(round(v, 2) * 100) % 25 == 0:
                print(f"   {v:12.2f} {a:14.2f}")
        # refine
        lo, hi = best - 0.05, best + 0.05
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if (ang(lo)[0] - 90) * (ang(mid)[0] - 90) <= 0: hi = mid
            else: lo = mid
        vbest = 0.5 * (lo + hi)
        a, dA, dB, db, cA, cB = ang(vbest)
        print(f"   -> 90 deg at v_D = {vbest:.3f} cm/us   (opening {a:.2f})")
        beamang = math.degrees(math.acos(np.clip(abs(db[2]), -1, 1)))
        print(f"      beam direction  = ({db[0]:+.3f},{db[1]:+.3f},{db[2]:+.3f})  -> {beamang:.2f} deg from z")
        print(f"      arm A direction = ({dA[0]:+.3f},{dA[1]:+.3f},{dA[2]:+.3f})")
        print(f"      arm B direction = ({dB[0]:+.3f},{dB[1]:+.3f},{dB[2]:+.3f})")
        angA = math.degrees(math.acos(np.clip(np.dot(dA, db), -1, 1)))
        angB = math.degrees(math.acos(np.clip(np.dot(dB, db), -1, 1)))
        print(f"      angles to beam  : A = {angA:.2f} deg,  B = {angB:.2f} deg,  sum = {angA+angB:.2f}")


if __name__ == "__main__":
    main()
