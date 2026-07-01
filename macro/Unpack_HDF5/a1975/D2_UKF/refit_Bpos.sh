#!/bin/bash
# Re-fit the 6 deuterium runs with B=+2.85 (sign test) reusing existing recos, 4 cores,
# output suffix _Bpos so the B=-2.85 fits are kept. Then ex_dp -> 17C Ex with the +2.85 fits.
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
OUT=/mnt/f/a1975/reco_d2/
LOGD=/tmp/claude-1000/-home-yassid/b5b91330-d765-4197-9f1c-6da8c3d5a215/scratchpad
echo "=== refit B=+2.85 (parallel 4 cores) $(date +%H:%M) ==="
for run in run_0300 run_0301 run_0302 run_0303 run_0304 run_0305; do
  root -l -b -q "fitGenfitter_a1975_deuterium.C(\"${run}_multifit\", 4000, \"$OUT\", \"_Bpos\", \"$OUT\", 2.85)" > $LOGD/fitBpos_${run}.log 2>&1 &
  while [ "$(jobs -r | wc -l)" -ge 4 ]; do sleep 3; done
done
wait
echo "=== ex_dp with the B=+2.85 fits $(date +%H:%M) ==="
RUNS="run_0300_multifit,run_0301_multifit,run_0302_multifit,run_0303_multifit,run_0304_multifit,run_0305_multifit"
root -l -b -q "ex_dp_a1975.C(\"$RUNS\",\"$OUT\",\"_Bpos\")" 2>&1 | grep -iE "candidates|17C Ex:|saved"
echo "REFIT_BPOS_DONE $(date +%H:%M)"
