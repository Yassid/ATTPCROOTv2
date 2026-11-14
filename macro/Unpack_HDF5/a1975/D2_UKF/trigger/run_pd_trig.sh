#!/bin/bash
# (p,d) trigger check, inputs: mesh features of the H2 runs + (p,d) track dump with entry index.
source ~/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
T=~/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/trigger; cd $T
seq 106 189 | xargs -P 4 -I{} sh -c 'n=$(printf %04d {}); root -l -b -q "dump_pd_tracks.C(\"run_$n\",\"/mnt/f/a1975/sim_rebuild/trig/pd/tracks_$n.txt\")" 2>&1 | grep done' > /mnt/f/a1975/sim_rebuild/trig/pd/dump.log
echo DUMP_DONE
seq 106 189 | xargs -P 6 -I{} ~/Spyral/venv/bin/python -I mesh_features_h2.py {} > /mnt/f/a1975/sim_rebuild/trig/mesh_h2/features.log 2>&1
echo FEATURES_DONE
