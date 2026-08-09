#!/usr/bin/env bash
# Reconstruct + fit the a1975 elastic simulation, so that reconstructed Ex can be plotted against
# RECONSTRUCTED vertex z -- the same quantities the data measurement uses. The truth-level slope
# is 0.487 MeV across the chamber and the data shows 0.05; this asks whether the reconstruction
# itself destroys it.
set -eo pipefail
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
HERE=$REPO/macro/Simulation/ATTPC/16C_pd
OUT=/mnt/f/a1975_C16_pp_chain; mkdir -p $OUT
set +u; source $REPO/build/config.sh >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
cd $HERE
echo "[$(date +%H:%M:%S)] reco"
root -b -q -l "run_reco_C16pd.C(\"/mnt/f/a1975_C16_pp_sim_el300.root\",\"$OUT/el300_reco.root\",\"ATTPC.a1975_C16_sim.par\",20,20,8,0,20.0,\"mover\")" > $OUT/reco.log 2>&1
grep -q "sim reco done" $OUT/reco.log || { echo "RECO_FAILED"; exit 1; }
echo "[$(date +%H:%M:%S)] fit (proton hypothesis, same config as the pp production)"
cd $REPO/macro/Unpack_HDF5/a1975/UKF
root -b -q -l "pipeline/fitGenfitter_a1975.C(\"el300\",-1,\"$OUT/\",\"\",\"$OUT/\",-2.85,2,5,\"proton_H2_catima.txt\",\"\",4.0,kFALSE,2212,1.00727647,1,kFALSE,\"ATTPC_H300torr_RT\")" > $OUT/fit.log 2>&1
[ -s "$OUT/el300_genfitter.root" ] && echo "[$(date +%H:%M:%S)] DONE" || echo "FIT_FAILED"
