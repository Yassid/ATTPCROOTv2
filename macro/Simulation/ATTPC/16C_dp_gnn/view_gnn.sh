#!/bin/bash
# Event-by-event viewer: HDBSCAN-trained GNN clustering (left) vs Spyral (right), matched events.
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC/16C_dp_gnn
~/gnn_env/bin/python build_gnn_viewer.py
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh 2>/dev/null
root -l 'viewer/view_events.C("data/cmp_gnn_0305.csv","data/cmp_spyral_0305.csv")'
