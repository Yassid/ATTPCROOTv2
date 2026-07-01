#!/bin/bash
# Working point: thr40 + AtPSAMultiFit + HDBSCAN(mcs20,ms8). Regenerate all 6 runs -> dump labels ->
# build parquet -> retrain GNN on HDBSCAN labels.
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh 2>/dev/null
D2=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
GNN=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC/16C_dp_gnn
cd $D2
echo "===== HDBSCAN PIPELINE START $(date) ====="
for r in 0300 0301 0302 0303 0304 0305; do
  echo "##### reco run_$r (thr40 + HDBSCAN) #####"
  root -l -b -q "unpackReco_multifit.C(\"run_$r\", 0, false, \"/mnt/f/a1975/reco_d2/\", \"/home/yassid/spyral_d2/h5/\", false, true, \"multifit\", 0, 40, \"hdbscan\", 20, 8)" 2>&1 | grep -iE "Done ->|Real time" | tail -2
done
echo "##### dump HDBSCAN labels (all 6) #####"
root -l -b -q dump_triplclust_all.C 2>&1 | grep -iE "run_0|wrote"
cd $GNN
echo "##### build parquet #####"
~/gnn_env/bin/python sup/build_attpc_parquet.py
echo "##### retrain GNN on HDBSCAN labels #####"
~/gnn_env/bin/python -u sup/train_attpc.py
echo "===== HDBSCAN PIPELINE DONE $(date) ====="
