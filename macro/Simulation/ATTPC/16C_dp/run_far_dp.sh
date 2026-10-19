#!/bin/bash
# 16C(d,p)17C g.s. with maxELoss 11.5 MeV so the generated vertices reach the far end (vz 900-1000), where the data
# has its forward-cm bound protons and the 10 MeV production (gs_s4001) has none. 2026-10-12.
# Per seed: generate (CM 2-50, Ebeam 184.25, B -28.5 kG) -> run_reco_dp.sh (gain 75000, data thr/mcs/ms) ->
# D2_UKF/census_arc_dp.C (Spyral + hit-helix per pattern track). No genfit (7 h / 20k events; Spyral suffices here).
#   ./run_far_dp.sh SEED NEV
SEED=$1; NEV=${2:-5000}; OUT=/mnt/f/a1975_C16_dp_sim_far; MAXEL=11.5
HERE=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC/16C_dp
set +u; source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export ROOT_INCLUDE_PATH="/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:$HOME/fair_install/FairRootInstall/include"
J=far_s${SEED}; SIM=$OUT/${J}_sim.root
cd $HERE
root -b -l -q "C16_dp_sim.C($NEV,2.0,50.0,\"TGeant4\",-28.5,\"$SIM\",0.0,$SEED,184.25,$MAXEL)" > $OUT/${J}_sim.log 2>&1 || true
[ -s $SIM ] || { echo "$J SIM FAILED"; exit 1; }
./run_reco_dp.sh $SIM 75000 || true
RECO=$OUT/${J}_reco_g75000.root
[ -s $RECO ] || { echo "$J RECO FAILED"; exit 1; }
rm -f $OUT/${J}_sim.root.tmp
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
root -l -b -q "census_arc_dp.C(\"$RECO\",\"$OUT/${J}_census.root\",1,$SEED)" > $OUT/${J}_census.log 2>&1 || true
echo "$J DONE"
