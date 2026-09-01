#!/usr/bin/env python3
"""Determine TBEntrance from the beam's closest approach to the detector axis.

TBEntrance fixes the longitudinal origin, and in thick-target inverse kinematics the vertex
position IS the beam energy -- so this number sets the energy scale of the excitation
function.

RESULT: TBEntrance = 320 +- 4, i.e. the value already in the parameter file is correct.

    run 100, magnet OFF        tb* = 320.0 +- 4.2   miss 12.5 mm
    run 128, corrected         tb* = 318.5 +- 10.5  miss 16.3 mm
    run 128, UNCORRECTED       tb* =  99.7 +- 25.1  miss 81.1 mm    <- control

The apparent "157 mm offset" that motivated this was a METHOD ARTEFACT. Two earlier
estimates (846 and 1157 mm) both extrapolated a closest-approach point PER EVENT from a
short lever arm, so a small slope error moved the crossing far; the per-event distribution
has an IQR of ~400 mm, and a median of fragile numbers is still fragile.

The control line is also the sharpest demonstration of the drift correction so far: without
it the beam's closest approach lands 220 tb (~790 mm) away from where the magnet-off run
puts it, and misses the axis by 81 mm. With it, 319 +- 10 against 320 +- 4 -- agreement to
one time bucket, with nothing tuned to achieve that.

METHOD: fit each event's beam line, then average the LINE PARAMETERS and evaluate the
closest approach once.

Do NOT instead average charge-weighted x,y in bins of tb and fit that. It looks more robust
and is badly biased: at a given tb the events do not all contain the same track, so the mean
is diluted toward the pad-plane centre and the fitted slope shrinks. That version returned a
beam polar angle of 3.72 deg where every other method gives ~6.4 -- classic regression
dilution -- and an implied shift of +580 mm.

CAVEAT: the polar angle from THIS estimator runs low (6.10 deg on run 100, 4.67 on run 128)
because taking a median over events that include some contaminated fits shrinks the slope.
Quote angles from the clustered fit in verify_production.py instead. The closest-approach
COMPARISON between runs is unaffected, since the same shrinkage applies to all three.

Run it on a MAGNET-OFF run for the reference value: with B = 0 there is no E x B shear, so
the measured trajectory is the beam itself and no drift correction has to be trusted.

usage: find_tbentrance.py [hits.csv] [--corrected]
"""
import sys, csv, math
from collections import defaultdict
import numpy as np

PATH = sys.argv[1] if len(sys.argv) > 1 else "/home/yassid/dec2014_calib/hits_run100.csv"
USE_COR = "--corrected" in sys.argv

TBTIME_US = 0.16          # 160 ns per time bucket
VD = 2.251                # cm/us
K = TBTIME_US * VD * 10   # mm per time bucket
ZPAD = 1000.0             # the z assigned to tb = TBEntrance by CalculateZGeo
TBENT_CURRENT = 320

NMIN, NMAX = 30, 120      # beam-like events: a reaction adds arms that drag the mean
BINW = 5


def load(path):
    ev = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh):
            x, y = (r["xc"], r["yc"]) if USE_COR else (r["x"], r["y"])
            ev[int(r["event"])].append((float(r["tb"]), float(x), float(y), float(r["q"])))
    return ev


ev = load(PATH)
rows = []
for h in ev.values():
    if NMIN <= len(h) <= NMAX:
        rows.extend(h)
a = np.array(rows)
if len(a) == 0:
    sys.exit("no beam-like events found")
tb, x, y, q = a[:, 0], a[:, 1], a[:, 2], a[:, 3]

print(f"{PATH}   {'corrected' if USE_COR else 'raw'} positions")
print(f"  {len(ev)} events, {len(a)} hits in beam-like events")
print(f"  pad-plane extent: x [{x.min():.0f}, {x.max():.0f}]  y [{y.min():.0f}, {y.max():.0f}] mm "
      f"-- the detector axis is x = y = 0")

# --- mean beam LINE ------------------------------------------------------------------
# Fit each event's beam line first, then average the LINE PARAMETERS -- not the positions.
#
# Averaging charge-weighted mean x,y in bins of tb and fitting that looked more robust and
# is badly biased: at a given tb the events do not all contain the same track, so the mean
# is diluted toward the pad-plane centre and the fitted slope shrinks. It returned a beam
# polar angle of 3.72 deg where every other method gives 6.4 deg -- classic regression
# dilution. Fitting per event is unbiased for the slope; averaging the parameters afterwards
# keeps the non-linear closest-approach step to a single evaluation, which was the reason
# for wanting an ensemble in the first place.
def fit_events(ev):
    P = []
    for h in ev.values():
        if not (NMIN <= len(h) <= NMAX):
            continue
        b = np.array(h)
        t, u, v, w = b[:, 0], b[:, 1], b[:, 2], b[:, 3]
        if np.ptp(t) < 50:               # need a lever arm in tb
            continue
        px = np.polyfit(t, u, 1, w=w)
        py = np.polyfit(t, v, 1, w=w)
        P.append((px[0], px[1], py[0], py[1]))
    return np.array(P)


def closest_approach(P):
    """Median line parameters -> single closest-approach evaluation."""
    bx, ax_, by, ay_ = np.median(P, axis=0)
    den = bx * bx + by * by
    t = -(bx * ax_ + by * ay_) / den
    miss = math.hypot(bx * t + ax_, by * t + ay_)
    return t, miss, bx, by


P = fit_events(ev)
print(f"  {len(P)} events with a usable beam line")
tstar, miss, bx, by = closest_approach(P)
print()
print(f"  dx/dtb = {bx:+.4f} mm/tb   dy/dtb = {by:+.4f} mm/tb")
print(f"  beam polar angle = {math.degrees(math.atan2(math.hypot(bx, by), K)):.2f} deg"
      f"   (sanity: must be ~6.4, not ~3.7)")
print()
print(f"  closest approach to the detector axis at tb = {tstar:.1f}")
print(f"  miss distance there                        = {miss:.1f} mm")
print()
print(f"  current TBEntrance = {TBENT_CURRENT}")
print(f"  implied shift      = {tstar - TBENT_CURRENT:+.1f} tb = "
      f"{(tstar - TBENT_CURRENT) * K:+.0f} mm in z")

# bootstrap over events
rng = np.random.default_rng(12345)
ts = []
for _ in range(400):
    idx = rng.choice(len(P), len(P), replace=True)
    t, _m, _a, _b = closest_approach(P[idx])
    ts.append(t)
ts = np.array(ts)
print()
print(f"  bootstrap (400 resamples): tb* = {ts.mean():.1f} +- {ts.std():.1f} "
      f"({ts.std()*K:.0f} mm in z)")
print(f"  => TBEntrance = {ts.mean():.0f} +- {ts.std():.0f}")
