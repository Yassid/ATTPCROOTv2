#!/bin/bash
# Full pipeline on the 6 deuterium runs, 4 cores: reco (thr40, HDBSCAN mover, full cloud)
# -> genfit fit (proton hyp, backward seed-fix). Each run's reco then fit runs as one job;
# up to 4 jobs concurrently. Outputs reco_d2/<run>_multifit_{reco,genfitter_p}.root.
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
OUT=/mnt/f/a1975/reco_d2/
H5=/home/yassid/spyral_d2/h5/
LOGD=/tmp/claude-1000/-home-yassid/b5b91330-d765-4197-9f1c-6da8c3d5a215/scratchpad
N=4000
recofit(){
  local run=$1
  root -l -b -q "unpackReco_multifit.C(\"$run\", $N, false, \"$OUT\", \"$H5\", false, true, \"multifit\", 0, 40, \"hdbscan\", 0, 3)" > $LOGD/log_${run}_reco.txt 2>&1
  root -l -b -q "fitGenfitter_a1975_deuterium.C(\"${run}_multifit\", $N, \"$OUT\", \"\", \"$OUT\", -2.85)" > $LOGD/log_${run}_fit.txt 2>&1
  echo "$(date +%H:%M) DONE $run"
}
for run in run_0300 run_0301 run_0302 run_0303 run_0304 run_0305; do
  recofit "$run" &
  while [ "$(jobs -r | wc -l)" -ge 4 ]; do sleep 5; done
done
wait
echo "ALL_RECOFIT_DONE $(date +%H:%M)"
