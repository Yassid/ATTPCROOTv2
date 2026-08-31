import os
#!/usr/bin/env python3
"""alpha+alpha elastic kinematics -- an acceptance-independent closure test.

Acceptance changes how densely events populate a kinematic locus; it does not move the
locus. So the elastic relations can be tested even though the efficiency is unknown:

    E1 = E_beam cos^2(theta1)      E2 = E_beam sin^2(theta1)
    theta1 + theta2 = 90 deg
    E1 + E2 = E_beam               (energy conservation)

Each quantity is measured a different way, which is what makes this a closure test rather
than a tautology:

  * E_beam  from the VERTEX POSITION via CATIMA energy loss
  * theta   from the TRACK DIRECTIONS at the vertex
  * E1, E2  from the TRACK LENGTHS via CATIMA range

Nothing here is shared with the calibration that produced the positions, so agreement is
a real check on the whole chain.
"""
import math, csv, sys, warnings, collections
import numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_tracks import cluster, fit_line, vertex_from_tracks, closest_approach

K = 0.16 * 2.251 * 10
TB0 = 42.3
ZW = 1000.0
PAD_R = 250.0        # pad-plane radius, mm -- used to reject tracks that leave the volume


def load_table(path="/home/yassid/dec2014_calib/eloss_alpha_heco2.txt"):
    d, e = [], []
    for line in open(path):
        if line.startswith("#"):
            continue
        a, b = line.split()
        d.append(float(a)); e.append(float(b))
    return np.array(d), np.array(e)


D, E = load_table()
RMAX = float(D[np.argmin(np.abs(E))])          # full range of a 7.8 MeV alpha


def energy_after(dist):
    """Beam energy after traversing `dist` mm of gas."""
    return float(np.interp(dist, D, E))


def range_of(energy):
    """Range of an alpha of the given energy: R(E) = Rmax - d(E)."""
    d_at = float(np.interp(energy, E[::-1], D[::-1]))
    return RMAX - d_at


def energy_from_range(r):
    """Invert R(E). Monotonic, so a direct interpolation on the tabulated curve."""
    grid = np.arange(0.05, 7.81, 0.05)
    rr = np.array([range_of(g) for g in grid])
    return float(np.interp(r, rr, grid))


def arc_length(P, q):
    """Path length along a track: order hits by projection on the fitted axis and sum."""
    c, d = fit_line(P, q)
    s = np.argsort(P @ d)
    Q = P[s]
    return float(np.sum(np.linalg.norm(np.diff(Q, axis=0), axis=1)))


def analyse(hits):
    A = np.array(hits)
    z = (A[:, 0] - TB0) * K
    P = np.column_stack([A[:, 1], A[:, 2], z])
    q = A[:, 3]
    lab, _, _ = cluster(P, q)
    segs = []
    for L in sorted(set(lab) - {-1}):
        m = lab == L
        if m.sum() < 8:
            continue
        c, d = fit_line(P[m], q[m])
        segs.append(dict(m=m, c=c, d=d, n=int(m.sum()),
                         L=float(np.ptp(P[m] @ d)), qm=float(q[m].mean())))
    if len(segs) < 3:
        return None
    cand = [s for s in segs if abs(s["d"][2]) > 0.85]
    if not cand:
        return None
    beam = min(cand, key=lambda s: s["qm"])
    arms = sorted([s for s in segs if s is not beam], key=lambda s: -s["L"])[:2]
    if len(arms) < 2:
        return None
    v = vertex_from_tracks([(a["c"], a["d"]) for a in arms] + [(beam["c"], beam["d"])])
    if not (0 < v[2] < ZW):
        return None
    if closest_approach(arms[0]["c"], arms[0]["d"], arms[1]["c"], arms[1]["d"])[1] > 30:
        return None

    db = beam["d"]
    if db[2] < 0:
        db = -db                      # beam points from the window toward the pad plane
    ebeam = energy_after(ZW - v[2])

    # Elastic selection. Without it the sample is dominated by the two failure
    # topologies: a single track split longitudinally (opening ~180 deg, which shows up
    # as theta1+theta2 = 180) and two parallel pile-up beams (~0 deg).
    o1 = arms[0]["d"] if np.dot(arms[0]["d"], arms[0]["c"] - v) > 0 else -arms[0]["d"]
    o2 = arms[1]["d"] if np.dot(arms[1]["d"], arms[1]["c"] - v) > 0 else -arms[1]["d"]
    opening = math.degrees(math.acos(np.clip(np.dot(o1, o2), -1, 1)))
    if not (55.0 < opening < 125.0):
        return None

    out = []
    for a in arms:
        Q = P[a["m"]]
        # require containment: the track must stop inside, not leave through a boundary
        rad = np.hypot(Q[:, 0], Q[:, 1]).max()
        if rad > PAD_R or Q[:, 2].min() < 20 or Q[:, 2].max() > ZW - 20:
            return None
        d = a["d"] if np.dot(a["d"], a["c"] - v) > 0 else -a["d"]
        th = math.degrees(math.acos(np.clip(np.dot(d, -db), -1, 1)))
        # angle to the beam direction of travel
        th = 180.0 - th if th > 90 else th
        rng = arc_length(Q, q[a["m"]])
        out.append((th, rng, energy_from_range(rng)))
    (t1, r1, e1), (t2, r2, e2) = out
    return dict(z=v[2], ebeam=ebeam, t1=t1, t2=t2, e1=e1, e2=e2, r1=r1, r2=r2)


def main():
    runs = [int(a) for a in sys.argv[1:]] or [128, 130, 134, 137]
    res = []
    for run in runs:
        ev = collections.defaultdict(list)
        for r in csv.DictReader(open(f"/home/yassid/dec2014_calib/prod_{run}.csv")):
            ev[int(r["event"])].append((int(r["tb"]), float(r["x"]), float(r["y"]), float(r["q"])))
        for hits in ev.values():
            try:
                a = analyse(hits)
            except Exception:
                a = None
            if a:
                res.append(a)
    print(f"{len(res)} contained elastic candidates from runs {runs}")
    np.save("/home/yassid/dec2014_calib/kin.npy",
            np.array([[r["z"], r["ebeam"], r["t1"], r["t2"], r["e1"], r["e2"]] for r in res]))


if __name__ == "__main__":
    main()
