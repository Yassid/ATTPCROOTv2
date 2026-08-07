#!/usr/bin/env python3
"""Look at a single Dec2014 alpha event, raw vs Langevin-corrected.

The stored AtHit position is the *uncorrected* one: pad (x0,y0) plus
    z = ZPadPlane - (TBEntrance - tb) * dzPerTB
with no tilt and no Lorentz correction (the block that would apply it is
commented out in AtPSASimple2.cxx). Here we redo the transform from the raw
(x0, y0, tb) using the Langevin drift vector, following eqs. (4)-(10) of the
AT-TPC commissioning paper (Bradt et al., NIM A 875 (2017) 65-79):

    omega*tau = (B/E) v_D
    v_x = v_D/(1+wt^2) * wt * sin(theta_t)
    v_y = v_D/(1+wt^2) * wt^2 * sin(theta_t) cos(theta_t)
    v_z = v_D/(1+wt^2) * (1 + wt^2 cos^2(theta_t))

    x = x0 - v_x * T/f ,  y = y0 - v_y * T/f ,  z = v_z * T/f

with T the drift time in time buckets measured from the pad plane.
"""
import csv, math, collections, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- detector / electronics constants (ATTPC.alpha_150torr.par) ---------------
B_FIELD = 0.5691       # T
E_FIELD = 5000.0       # V/m
TILT    = 7.0          # deg
V_D     = 2.25         # cm/us   <-- the number we are trying to calibrate
TB_NS   = 160.0        # ns per time bucket (SamplingRate 6 -> 160 ns)
TB_ENTRANCE = 280      # tb of the z reference used by CalculateZGeo
Z_PADPLANE  = 1000.0   # mm


def drift_vector(v_d=V_D, b=B_FIELD, e=E_FIELD, tilt_deg=TILT):
    """Langevin drift velocity components, cm/us. v_d in cm/us."""
    t = math.radians(tilt_deg)
    wt = (b / e) * (v_d * 1e4)          # v_d cm/us -> m/s
    front = v_d / (1.0 + wt * wt)
    vx = front * wt * math.sin(t)
    vy = front * wt * wt * math.sin(t) * math.cos(t)
    vz = front * (1.0 + wt * wt * math.cos(t) ** 2)
    return vx, vy, vz, wt


def reconstruct(x0, y0, tb, v_d=V_D, tb_ref=TB_ENTRANCE, apply_lorentz=True):
    """Raw pad coords + time bucket -> calibrated position (mm)."""
    vx, vy, vz, _ = drift_vector(v_d)
    # drift time in us, measured from the pad-plane reference
    T = (tb_ref - np.asarray(tb, dtype=float)) * TB_NS * 1e-3
    if apply_lorentz:
        x = x0 - vx * T * 10.0      # cm/us * us -> cm -> mm
        y = y0 - vy * T * 10.0
        z = Z_PADPLANE - vz * T * 10.0
    else:
        x = np.asarray(x0, dtype=float)
        y = np.asarray(y0, dtype=float)
        z = Z_PADPLANE - v_d * T * 10.0
    return x, y, z


def load(path, event):
    hits = []
    for r in csv.DictReader(open(path)):
        if int(r["event"]) == event:
            hits.append((int(r["tb"]), float(r["x"]), float(r["y"]), float(r["q"])))
    tb = np.array([h[0] for h in hits])
    x0 = np.array([h[1] for h in hits])
    y0 = np.array([h[2] for h in hits])
    q = np.array([h[3] for h in hits])
    return tb, x0, y0, q


def main():
    event = int(sys.argv[1]) if len(sys.argv) > 1 else 258
    tb, x0, y0, q = load("/home/yassid/dec2014_calib/events_run0128.csv", event)
    print(f"event {event}: {len(tb)} hits")

    vx, vy, vz, wt = drift_vector()
    print(f"  omega*tau={wt:.3f}  v=({vx:.4f},{vy:.4f},{vz:.4f}) cm/us")

    xr, yr, zr = reconstruct(x0, y0, tb, apply_lorentz=False)
    xc, yc, zc = reconstruct(x0, y0, tb, apply_lorentz=True)
    print(f"  max Lorentz shift: dx={np.abs(xc-xr).max():.1f} mm  dy={np.abs(yc-yr).max():.1f} mm")

    s = 12 + 60 * (q / q.max())
    fig, ax = plt.subplots(2, 3, figsize=(16, 9))
    fig.suptitle(f"run_0128  event {event}   ({len(tb)} hits)   "
                 f"top: no Lorentz    bottom: Langevin-corrected "
                 f"(v_D={V_D} cm/us, tilt={TILT} deg)", fontsize=12)

    for row, (X, Y, Z, tag) in enumerate([(xr, yr, zr, "raw"), (xc, yc, zc, "corrected")]):
        for col, (A, B_, la, lb) in enumerate([
                (Z, X, "z [mm]", "x [mm]"),
                (Z, Y, "z [mm]", "y [mm]"),
                (X, Y, "x [mm]", "y [mm]")]):
            a = ax[row][col]
            sc = a.scatter(A, B_, c=q, s=s, cmap="viridis", norm=matplotlib.colors.LogNorm())
            a.set_xlabel(la); a.set_ylabel(lb)
            a.set_title(f"{tag}: {lb.split()[0]} vs {la.split()[0]}")
            a.grid(alpha=.3)
            if col == 2:
                a.set_aspect("equal")
    fig.colorbar(sc, ax=ax, label="charge", fraction=.02)
    out = f"/home/yassid/dec2014_calib/event_{event}.png"
    fig.savefig(out, dpi=110, bbox_inches="tight")
    print("wrote", out)


if __name__ == "__main__":
    main()
