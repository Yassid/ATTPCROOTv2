#!/bin/bash
# Elastic-deuteron embedding into real D2 events (embed_dd.C), 2026-11-10. Read run.log for EMBED_ALLDONE.
source ~/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export ROOT_INCLUDE_PATH=$HOME/fair_install/ATTPCROOTv2-OpenKF/build/include:$HOME/fair_install/FairRootInstall/include
M=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/trigger/dd
E=/mnt/f/a1975/caches/dd_r280/embed; SR=/mnt/f/a1975/sim_rebuild/r280/dd_el
cd $E
for rs in 16:6001 31:6002 79:6001 102:6002; do r=${rs%:*}; s=${rs#*:}
  root -l -b -q "$M/embed_dd.C($r,2,\"$E/e_r${r}_s${s}.root\",\"$SR/el_s${s}_reco.root\")" </dev/null > $E/e_r${r}_s${s}.log 2>&1 &
done
wait; grep -H "done ->" $E/e_*.log; echo EMBED_ALLDONE
