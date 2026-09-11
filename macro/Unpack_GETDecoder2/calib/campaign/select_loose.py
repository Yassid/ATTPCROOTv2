#!/usr/bin/env python3
"""Candidate selection with the arm cuts as PARAMETERS.

calibrate_run.select() hardcodes MIN_ARM_LEN = 60 mm and MIN_ARM_HITS = 12, which puts a hard
floor at 2.05 MeV per alpha -- i.e. E_cm >~ 2.05 MeV, right on top of the 8Be 2+ at 2.94.
The paper's own detector limit is ~200 keV / ~20 mm, so most of that is thrown away by choice,
not by the detector.

  usage: select_loose.py <slim.csv> <out.pkl> [min_arm_len] [min_arm_hits]
"""
import sys, os, pickle, warnings
import numpy as np
warnings.filterwarnings("ignore")
sys.path.insert(0, "/home/yassid/fair_install/ATTPCROOTv2_fr19port/macro/Unpack_GETDecoder2/calib")
from cluster_tracks import cluster, fit_line, closest_approach
from calibrate_run import load_events, TB_NS, TB0

MIN_LEN  = float(sys.argv[3]) if len(sys.argv) > 3 else 25.0
MIN_HITS = int(sys.argv[4])   if len(sys.argv) > 4 else 6
MAX_VTX, BEAM_ALIGN = 40.0, 0.85

def select(hits, v_ref=2.25):
    tb, x0, y0, q = hits[:, 0], hits[:, 1], hits[:, 2], hits[:, 3]
    P = np.column_stack([x0, y0, (tb - TB0) * TB_NS * 1e-3 * v_ref * 10.0])
    lab, _, _ = cluster(P, q)
    segs = []
    for L in sorted(set(lab) - {-1}):
        m = lab == L
        if m.sum() < 4:
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
    if len(arms) < 2 or any(a["n"] < MIN_HITS or a["length"] < MIN_LEN for a in arms):
        return None
    if closest_approach(arms[0]["c"], arms[0]["d"], arms[1]["c"], arms[1]["d"])[1] > MAX_VTX:
        return None
    return dict(tb=tb, x=x0, y=y0, q=q, beam=beam["mask"], a1=arms[0]["mask"], a2=arms[1]["mask"])

csv_path, out_path = sys.argv[1], sys.argv[2]
evs, groups = load_events(csv_path)
cands = []
for h in groups:
    try:
        r = select(h)
    except Exception:
        r = None
    if r:
        cands.append(r)
pickle.dump({"nevents": len(evs), "cands": cands}, open(out_path, "wb"))
print(f"{os.path.basename(csv_path)}: events={len(evs)} candidates={len(cands)} "
      f"(min_len={MIN_LEN:.0f} mm, min_hits={MIN_HITS})")
