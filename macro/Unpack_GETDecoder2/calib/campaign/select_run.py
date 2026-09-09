#!/usr/bin/env python3
"""Select alpha+alpha elastic candidates from one run's slim hit CSV.

Reuses calibrate_run.select so the campaign-wide numbers come from exactly the same
selection that produced the single-run calibration (v_D = 2.251, tilt = 6.47 on run 100).
Candidates keep their hits, so the geometry can be recomputed for any (v_D, E, tilt)
afterwards without re-clustering.

  usage: select_run.py <slim.csv> <out.pkl>
"""
import sys, os, pickle, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, "/home/yassid/fair_install/ATTPCROOTv2_fr19port/macro/Unpack_GETDecoder2/calib")
from calibrate_run import load_events, select

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
with open(out_path, "wb") as fh:
    pickle.dump({"nevents": len(evs), "cands": cands}, fh)
print(f"{os.path.basename(csv_path)}: events={len(evs)} candidates={len(cands)}")
