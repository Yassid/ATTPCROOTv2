#!/bin/bash
# Which drift velocity makes the two energy determinations consistent?
#
# Each event carries two independent E_cm: from the two arm RANGES, and from how deep the
# beam got.  The tuned simulation says what their difference should be once reconstruction
# bias is taken out (arms -0.089, vertex +0.481, so vertex minus arms = +0.570 MeV).  The
# data give -0.276, an 0.85 MeV disagreement.
#
# v_D is the natural handle: it stretches z, and the beam path is nearly along z while the
# arms sit near 45 degrees, so raising v_D lengthens the path much more than the arms.  A
# longer path leaves less residual beam energy, pushing the vertex estimate DOWN relative
# to the arms -- which is the direction of the discrepancy.
cd /home/yassid/dec2014_package
export DEC2014_DATA=/home/yassid/dec2014_retrack
echo "[$(date +%T)] v_D scan starting"
for VD in 2.100 2.142 2.180 2.237 2.278; do
  DEC2014_VD=$VD python3 - <<'PY'
import os, sys
sys.path.insert(0, '/home/yassid/dec2014_package')
import numpy as np, config
from lib import loader
d = loader.build(config, quiet=True)
CEIL = config.E_BEAM / 2.0
Ecm = (d["E1"] + d["E2"]) / 2.0
keep = (np.abs(d["thsum"] - 90) < config.SUM_WINDOW) & (Ecm <= CEIL)
from lib.eloss import RangeTable
tab = RangeTable(config.ELOSS, config.PRESSURE, config.P_TABLE)
vtx = tab.energy_after(d["path"][keep]) / 2.0
arm = Ecm[keep]
dif = np.median(vtx - arm)
# the simulation's expectation for the same quantity
SIM = 0.570
print(f"  v_D {config.V_DRIFT:.3f} : N={keep.sum():6d}  median(vtx-arms) = {dif:+.3f} MeV"
      f"   gap to simulation = {dif-SIM:+.3f}   median E_cm(arms) {np.median(arm):.2f}")
PY
done
echo "[$(date +%T)] v_D scan done"
