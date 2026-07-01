#!/bin/bash
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh 2>/dev/null
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
RC(){ root -l -b -q "recluster.C(\"$1\",\"$2\",\"$3\",$4,$5,$6,$7)" 2>&1 | grep -iE "wrote" | tail -1; }
for THR in t40 t20; do
  IN=/tmp/claude-1000/-home-yassid/b5b91330-d765-4197-9f1c-6da8c3d5a215/scratchpad/${THR}_run_0305_multifit_reco.root
  RC $IN /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC/16C_dp_gnn/data/sw_${THR}_tc_a03_t4.csv   tc 0.03 4  15 6
  RC $IN /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC/16C_dp_gnn/data/sw_${THR}_tc_a03_t8.csv   tc 0.03 8  15 6
  RC $IN /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC/16C_dp_gnn/data/sw_${THR}_tc_a03_t12.csv  tc 0.03 12 15 6
  RC $IN /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC/16C_dp_gnn/data/sw_${THR}_tc_a10_t8.csv   tc 0.10 8  15 6
  RC $IN /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC/16C_dp_gnn/data/sw_${THR}_tc_a20_t12.csv  tc 0.20 12 15 6
  RC $IN /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC/16C_dp_gnn/data/sw_${THR}_hdbscan.csv     hdbscan 0 0 15 6
done
echo SWEEP_DONE
