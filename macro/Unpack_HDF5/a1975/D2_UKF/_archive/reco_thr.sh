#!/bin/bash
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh 2>/dev/null
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
# thr40, prominence off, 1500 events -> scratchpad
root -l -b -q 'unpackReco_multifit.C("run_0305", 1500, false, "/tmp/claude-1000/-home-yassid/b5b91330-d765-4197-9f1c-6da8c3d5a215/scratchpad/t40_", "/home/yassid/spyral_d2/h5/", false, true, "multifit", 0, 40)'
# thr20 same config for comparison
root -l -b -q 'unpackReco_multifit.C("run_0305", 1500, false, "/tmp/claude-1000/-home-yassid/b5b91330-d765-4197-9f1c-6da8c3d5a215/scratchpad/t20_", "/home/yassid/spyral_d2/h5/", false, true, "multifit", 0, 20)'
echo RECO_THR_DONE
