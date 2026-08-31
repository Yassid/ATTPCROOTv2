import os
#!/usr/bin/env python3
"""Cluster every event of run_0128 and collect alpha+alpha elastic candidates.

Clustering is done once per event at a reference drift velocity; the cluster
assignment is then reused while v_D is scanned, so HDBSCAN runs 11k times rather
than 11k x N_scan times. That is safe over a modest v_D range because scaling z
stretches the tracks but does not reorder which hit belongs to which arm.

For each surviving event we store the hits of the three segments (beam + two
arms) so the geometry can be recomputed for any (v_D, E, tilt) afterwards.
"""
import csv, math, sys, pickle, warnings
import numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cluster_tracks import cluster, fit_line, vertex_from_tracks, closest_approach

V_REF = 2.25
TB_NS = 160.0
TB0 = 2.2

MIN_ARM_HITS = 12
MIN_ARM_LEN = 60.0      # mm
MAX_VTX_DIST = 40.0     # mm, closest approach of the two arms
BEAM_ALIGN = 0.85       # |dz| for the beam candidate


def load_all(path):
    ev, tb, x, y, q = [], [], [], [], []
    with open(path) as fh:
        for r in csv.reader(fh):
            if r[0] == "event":
                continue
            ev.append(int(r[0])); tb.append(int(r[1]))
            x.append(float(r[2])); y.append(float(r[3])); q.append(float(r[4]))
    ev = np.array(ev); order = np.argsort(ev, kind="stable")
    ev, tb = ev[order], np.array(tb)[order]
    x, y, q = np.array(x)[order], np.array(y)[order], np.array(q)[order]
    bounds = np.searchsorted(ev, np.unique(ev))
    return np.unique(ev), np.split(np.column_stack([tb, x, y, q]), bounds[1:])


def analyse(hits):
    tb, x0, y0, q = hits[:, 0], hits[:, 1], hits[:, 2], hits[:, 3]
    T = (tb - TB0) * TB_NS * 1e-3
    P = np.column_stack([x0, y0, V_REF * T * 10.0])
    lab, _, _ = cluster(P, q)
    segs = []
    for L in sorted(set(lab) - {-1}):
        m = lab == L
        if m.sum() < 6:
            continue
        c, d = fit_line(P[m], q[m])
        segs.append(dict(mask=m, c=c, d=d, n=int(m.sum()),
                         length=float(np.ptp(P[m] @ d)), qmean=float(q[m].mean())))
    if len(segs) < 3:
        return None
    # beam: well aligned with z and the lowest mean charge among those
    cand = [s for s in segs if abs(s["d"][2]) > BEAM_ALIGN]
    if not cand:
        return None
    beam = min(cand, key=lambda s: s["qmean"])
    arms = sorted([s for s in segs if s is not beam], key=lambda s: -s["length"])[:2]
    if len(arms) < 2:
        return None
    if any(a["n"] < MIN_ARM_HITS or a["length"] < MIN_ARM_LEN for a in arms):
        return None
    _, dist = closest_approach(arms[0]["c"], arms[0]["d"], arms[1]["c"], arms[1]["d"])
    if dist > MAX_VTX_DIST:
        return None
    return dict(tb=tb, x=x0, y=y0, q=q,
                beam=beam["mask"], a1=arms[0]["mask"], a2=arms[1]["mask"],
                vtx_dist=float(dist))


def main():
    evs, groups = load_all("/home/yassid/dec2014_calib/all_run0128.csv")
    print(f"loaded {len(evs)} events", flush=True)
    out, nfail = {}, 0
    for i, (e, h) in enumerate(zip(evs, groups)):
        if i % 500 == 0:
            print(f"  {i}/{len(evs)}  kept={len(out)}", flush=True)
        try:
            r = analyse(h)
        except Exception:
            r = None
        if r is None:
            nfail += 1
        else:
            out[int(e)] = r
    print(f"\nselected {len(out)} candidate events of {len(evs)} ({nfail} rejected)")
    with open("/home/yassid/dec2014_calib/candidates.pkl", "wb") as fh:
        pickle.dump(out, fh)
    print("wrote candidates.pkl")


if __name__ == "__main__":
    main()
