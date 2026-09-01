#!/usr/bin/env python3
"""Figures for the manual: the beam-frame view and the validation summary.

Figure 1 reproduces what the event display shows (raw / corrected / de-tilted hits plus the
two reference axes) but as a 2D projection with a common scale, because in the TEve view the
three point sets differ by only a few degrees and visually coincide.

Figure 2 collects the validation results in one place: the MC-truth direction test binned by
drift, and the rotate-to-beam-axis test on real data.

usage: make_figures.py [hits_run128.csv] [outdir]
"""
import sys, csv, math
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

SRC = sys.argv[1] if len(sys.argv) > 1 else "/home/yassid/dec2014_calib/hits_run128.csv"
OUT = sys.argv[2] if len(sys.argv) > 2 else "/home/yassid/dec2014_calib/plots"

# Beam axis MEASURED from run 100 (B=0). Not the parameter-file TiltAng/ThetaRot: rotating
# by the measured values leaves the corrected beam 0.99 deg off the detector axis, the
# parameter values 1.91 deg (manual 5.4.8).
TILT, AZIM = 5.44, -173.2


def rot_matrix(tilt, azim):
    t, p = math.radians(tilt), math.radians(azim)
    b = np.array([math.sin(t) * math.cos(p), math.sin(t) * math.sin(p), math.cos(t)])
    z = np.array([0.0, 0.0, 1.0])
    v, c = np.cross(b, z), float(b @ z)
    s = np.linalg.norm(v)
    K = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]])
    return np.eye(3) + K + K @ K * ((1 - c) / s**2), b


R, BDIR = rot_matrix(TILT, AZIM)

ev = defaultdict(list)
with open(SRC) as fh:
    for r in csv.DictReader(fh):
        ev[int(r["event"])].append((float(r["x"]), float(r["y"]), float(r["xc"]),
                                    float(r["yc"]), float(r["z"]), float(r["q"])))

# Pick a long, clean, beam-like event so the geometry is legible.
best, bestspan = None, 0
for k, h in ev.items():
    if not (40 <= len(h) <= 110):
        continue
    a = np.array(h)
    span = a[:, 4].max() - a[:, 4].min()
    if span > bestspan:
        best, bestspan = k, span
a = np.array(ev[best])
raw = np.column_stack([a[:, 0], a[:, 1], a[:, 4]])
cor = np.column_stack([a[:, 2], a[:, 3], a[:, 4]])
det = cor @ R.T

fig, ax = plt.subplots(1, 2, figsize=(15, 5.6))


def mean_line(P_of):
    """Mean slope and intercept over events -- the ensemble trajectory.

    A single event is NOT representative: per-event direction scatters by ~5 deg (beam
    divergence plus fit noise), so one event's de-tilted track is not flat even though the
    ensemble is. Drawing one event here would misrepresent the result.
    """
    S, I = [], []
    for h in ev.values():
        if not (30 <= len(h) <= 120):
            continue
        b = np.array(h)
        if b[:, 4].max() - b[:, 4].min() < 300:
            continue
        P = P_of(b)
        zc = P[:, 2] - P[:, 2].mean()
        den = (zc * zc).sum()
        if den <= 0:
            continue
        sx = (zc * (P[:, 0] - P[:, 0].mean())).sum() / den
        S.append(sx)
        I.append(P[:, 0].mean() - sx * P[:, 2].mean())
    return float(np.mean(S)), float(np.median(I))


a0 = ax[0]
zl = np.array([0.0, 1200.0])
cases = (
    (lambda b: np.column_stack([b[:, 0], b[:, 1], b[:, 4]]), "tab:red", "raw (uncorrected)"),
    (lambda b: np.column_stack([b[:, 2], b[:, 3], b[:, 4]]), "tab:blue", "Lorentz-corrected"),
    (lambda b: np.column_stack([b[:, 2], b[:, 3], b[:, 4]]) @ R.T, "tab:green",
     "corrected, then de-tilted"),
)
# faint hits of one event for context only
a0.plot(raw[:, 2], raw[:, 0], ".", ms=3, color="0.8", alpha=.5, label="hits of one event (context)")
x0 = None
for f, c, lab in cases:
    sx, ic = mean_line(f)
    if x0 is None:
        x0 = ic                      # anchor the reference axes at the same start point,
    a0.plot(zl, ic + sx * zl, "-", color=c, lw=2.4, label=lab)   # so slopes compare by eye
# Reference directions drawn FROM THE BEAM, not from the coordinate origin: only the slope
# is meaningful, and anchoring them elsewhere makes the comparison unreadable.
a0.plot(zl, x0 + 0 * zl, "--", color="0.4", lw=1.4, label="detector-axis direction")
a0.plot(zl, x0 + BDIR[0] / BDIR[2] * zl, "--", color="tab:orange", lw=1.4,
        label="beam direction (= B), measured")
a0.set_xlabel("z  (drift direction) [mm]")
a0.set_ylabel("x [mm]   (x projection only)")
a0.set_title("Mean beam trajectory over ~2200 events, x projection\n"
             "the correction puts the beam on B; de-tilting then puts it on the detector axis")
a0.legend(fontsize=8, loc="best")
a0.grid(alpha=.3)

# --- ensemble angles, the quantitative version of the left panel ----------------------
def ens(P_of):
    S = []
    for h in ev.values():
        if not (30 <= len(h) <= 120):
            continue
        b = np.array(h)
        if b[:, 4].max() - b[:, 4].min() < 300:
            continue
        P = P_of(b)
        zc = P[:, 2] - P[:, 2].mean()
        den = (zc * zc).sum()
        if den <= 0:
            continue
        S.append(((zc * (P[:, 0] - P[:, 0].mean())).sum() / den,
                  (zc * (P[:, 1] - P[:, 1].mean())).sum() / den))
    sx, sy = np.array(S).mean(axis=0)
    d = np.array([sx, sy, 1.0])
    return math.degrees(math.acos(abs(d[2] / np.linalg.norm(d))))


vals = [
    ("raw\nunrotated", ens(lambda b: np.column_stack([b[:, 0], b[:, 1], b[:, 4]]))),
    ("corrected\nunrotated", ens(lambda b: np.column_stack([b[:, 2], b[:, 3], b[:, 4]]))),
    ("corrected\nde-tilted", ens(lambda b: np.column_stack([b[:, 2], b[:, 3], b[:, 4]]) @ R.T)),
]
a1 = ax[1]
cols = ["tab:red", "tab:blue", "tab:green"]
a1.bar(range(3), [v for _, v in vals], color=cols, width=.55)
for i, (_, v) in enumerate(vals):
    a1.text(i, v + .08, f"{v:.2f}°", ha="center", fontsize=11, weight="bold")
a1.axhline(TILT, color="tab:orange", ls="--", lw=1.3)
a1.text(2.45, TILT + .07, "beam tilt\n(measured, B=0)", fontsize=8, color="tab:orange", ha="right")
a1.set_xticks(range(3))
a1.set_xticklabels([k for k, _ in vals])
a1.set_ylabel("ensemble beam angle from the detector axis [deg]\n(full 3D angle, not the x projection at left)")
a1.set_ylim(0, 6.4)
a1.set_title("Corrected and de-tilted, the beam lands on the axis\n"
             "0.99° residual — the event looks as though the detector were never tilted")
a1.grid(alpha=.3, axis="y")

fig.tight_layout()
fig.savefig(f"{OUT}/beam_frame.png", dpi=110, bbox_inches="tight")
print(f"wrote {OUT}/beam_frame.png  (event {best}, {len(a)} hits)")
for k, v in vals:
    print(f"  {k.replace(chr(10),' '):22s} {v:.2f} deg")
