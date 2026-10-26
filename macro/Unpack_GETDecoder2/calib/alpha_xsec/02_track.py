#!/usr/bin/env python3
"""Step 2 -- TRACKING.  Find the beam and the two scattered tracks in each event.

Reads the CSVs made by 01_dump.sh, clusters each event's hits into tracks, keeps the
events that look like a reaction (a beam-like segment plus two others), and writes one
pickle per run.

Plain HDBSCAN on (x,y,z) merges the two outgoing tracks at the vertex, which is exactly
where you least want it to -- the vertex is a real density maximum where both tracks
overlap.  What separates them there is not position but DIRECTION, so each hit is given a
local direction from a PCA over its neighbours and the clustering runs in
[position, orientation] space.  The vertex is then recovered as the intersection of the
fitted lines, not as a cluster of its own, which is what makes it well defined.

NO arm-length cut is applied.  The clustering's own min_cluster_size already imposes a
floor near 5 mm, and every length cut this analysis inherited turned out to be deleting
physics rather than background.

    python3 02_track.py            all runs
    python3 02_track.py 3          just the first 3
"""
import sys, os, glob, pickle, csv, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np
import config
from lib.clustering import cluster, fit_line, closest_approach

V_REF       = 2.25     # reference drift velocity used ONLY to cluster; the physics
                       # downstream recomputes z from tb with config.V_DRIFT
BEAM_ALIGN  = 0.85     # a beam-like segment points within ~32 deg of the drift axis
MAX_VTX_MM  = 40.0     # the two tracks must approach within this to count as one vertex
MIN_SEG_HIT = 4        # a segment needs at least this many hits to be fitted at all


def load_csv(path):
    """-> (event ids, list of per-event arrays [tb, x, y, q])"""
    ev, tb, x, y, q = [], [], [], [], []
    with open(path) as fh:
        for r in csv.DictReader(fh):
            ev.append(int(r["event"])); tb.append(float(r["tb"]))
            x.append(float(r["x"])); y.append(float(r["y"])); q.append(float(r["q"]))
    ev = np.asarray(ev)
    order = np.argsort(ev, kind="stable")
    ev = ev[order]
    data = np.column_stack([np.asarray(tb)[order], np.asarray(x)[order],
                            np.asarray(y)[order], np.asarray(q)[order]])
    uniq = np.unique(ev)
    return uniq, np.split(data, np.searchsorted(ev, uniq)[1:])


def find_tracks(hits):
    """One event -> dict of masks, or None if it is not a reaction candidate."""
    tb, x, y, q = hits[:, 0], hits[:, 1], hits[:, 2], hits[:, 3]
    z = (tb - config.T0_BUCKETS) * config.TB_NS * 1e-3 * V_REF * 10.0
    P = np.column_stack([x, y, z])
    lab, _, _ = cluster(P, q)
    segs = []
    for L in sorted(set(lab) - {-1}):
        m = lab == L
        if m.sum() < MIN_SEG_HIT:
            continue
        c, d = fit_line(P[m], q[m])
        segs.append(dict(mask=m, c=c, d=d, n=int(m.sum()),
                         length=float(np.ptp(P[m] @ d)), qmean=float(q[m].mean())))
    if len(segs) < 3:
        return None
    beamish = [s for s in segs if abs(s["d"][2]) > BEAM_ALIGN]
    if not beamish:
        return None
    beam = min(beamish, key=lambda s: s["qmean"])      # the beam ionises least per pad
    arms = sorted([s for s in segs if s is not beam], key=lambda s: -s["length"])[:2]
    if len(arms) < 2:
        return None
    if closest_approach(arms[0]["c"], arms[0]["d"], arms[1]["c"], arms[1]["d"])[1] > MAX_VTX_MM:
        return None
    return dict(tb=tb, x=x, y=y, q=q,
                beam=beam["mask"], a1=arms[0]["mask"], a2=arms[1]["mask"])


if __name__ == "__main__":
    nmax = int(sys.argv[1]) if len(sys.argv) > 1 else 0
    files = sorted(glob.glob(os.path.join(config.CSV_DIR, "hits_*.csv")))
    if not files:
        raise SystemExit(f"no CSVs in {config.CSV_DIR} -- run ./01_dump.sh first")
    if nmax:
        files = files[:nmax]
    os.makedirs(config.TRACKS_DIR, exist_ok=True)
    for path in files:
        run = os.path.basename(path).replace("hits_", "").replace(".csv", "")
        out = os.path.join(config.TRACKS_DIR, f"tracks_{run}.pkl")
        if os.path.exists(out):
            print(f"  skip {run} (already tracked)"); continue
        evs, groups = load_csv(path)
        found = []
        for h in groups:
            try:
                r = find_tracks(h)
            except Exception:
                r = None
            if r:
                found.append(r)
        pickle.dump({"nevents": len(evs), "cands": found}, open(out, "wb"))
        print(f"  run {run}: {len(evs)} events -> {len(found)} candidates "
              f"({100*len(found)/max(len(evs),1):.0f}%)")
    print(f"  tracks in {config.TRACKS_DIR}")
