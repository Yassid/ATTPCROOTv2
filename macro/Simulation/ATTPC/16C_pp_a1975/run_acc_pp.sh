#!/usr/bin/env bash
# Generate -> reco -> fit -> acceptance for the a1975 (p,p) elastic, at the correct 300 torr gas.
set -eo pipefail
NEV=${1:-12000}; SEED=${2:-1002}
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
HERE=$REPO/macro/Simulation/ATTPC/16C_pp_a1975
OUT=/mnt/f/a1975_C16_pp_acc; mkdir -p $OUT
set +u; source $REPO/build/config.sh >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
J=el_s$SEED
cd $HERE
echo "[$(date +%H:%M:%S)] generate ($NEV events)"
root -b -q -l "C16_pp_a1975_sim.C($NEV,2.,178.,\"TGeant4\",-28.5,\"$OUT/${J}_sim.root\",0.0,$SEED)" > $OUT/${J}_gen.log 2>&1
grep -q "Macro finished successfully" $OUT/${J}_gen.log || { echo GEN_FAILED; exit 1; }
echo "[$(date +%H:%M:%S)] reconstruct"
cd $REPO/macro/Simulation/ATTPC/16C_pd
root -b -q -l "run_reco_C16pd.C(\"$OUT/${J}_sim.root\",\"$OUT/${J}_reco.root\",\"ATTPC.a1975_C16_sim.par\",20,20,8,0,20.0,\"mover\")" > $OUT/${J}_reco.log 2>&1
grep -q "sim reco done" $OUT/${J}_reco.log || { echo RECO_FAILED; exit 1; }
echo "[$(date +%H:%M:%S)] fit (proton hypothesis)"
cd $REPO/macro/Unpack_HDF5/a1975/UKF
root -b -q -l "pipeline/fitGenfitter_a1975.C(\"$J\",-1,\"$OUT/\",\"\",\"$OUT/\",-2.85,2,5,\"proton_H2_catima.txt\",\"\",4.0,kFALSE,2212,1.00727647,1,kFALSE,\"ATTPC_H300torr_RT\")" > $OUT/${J}_fit.log 2>&1
[ -s "$OUT/${J}_genfitter.root" ] || { echo FIT_FAILED; exit 1; }
echo "[$(date +%H:%M:%S)] acceptance"
cd $HERE
root -b -q -l "acceptance_C16pp.C(\"$OUT/${J}_sim.root\",\"$OUT/${J}_genfitter.root\",\"el\",0.0,188.0,5.0,36,180.0,10.0,0.5,2.0,kFALSE)" > $OUT/${J}_acc.log 2>&1
grep "overall acceptance" $OUT/${J}_acc.log
echo "[$(date +%H:%M:%S)] ACCDONE"
