#!/bin/bash
# Event-by-event ATTPCROOT-vs-Spyral cluster viewer.  Usage: bash view_compare.sh [run] [offset]
#   run    = run number (default 305)
#   offset = first GET event # of the chunk (default 34675 for run_0305; from h5 evt<N>_data keys)
RUN=${1:-305}; OFF=${2:-34675}
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC/16C_dp_gnn
~/gnn_env/bin/python build_compare_viewer.py $RUN $OFF
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh 2>/dev/null
root -l "viewer/view_events.C(\"data/cmp_attpc_0${RUN}.csv\",\"data/cmp_spyral_0${RUN}.csv\")"
