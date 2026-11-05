#!/bin/bash
# Pulse model for every r280 g.s. sim event with a genfit row, and truth key.
source /mnt/f/a1975/sim_rebuild/r280/env.sh; cd /mnt/f/a1975/sim_rebuild/trig
root -l -b -q '/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/geom_vars_dp.C("/mnt/f/a1975/sim_rebuild/r280/p_gs_sim.root","/mnt/f/a1975/sim_rebuild/trig/truth_r280.root")' </dev/null 2>&1 | tail -1
root -l -b -q 'dump_keys.C' </dev/null 2>&1 | tail -1
awk '{print $1}' keys.txt | sort -n -u > qhits/ev_all.txt
root -l -b -q 'mesh_q_hits.C("/mnt/f/a1975/sim_rebuild/r280/p_gs_reco.root","qhits/ev_all.txt","qhits/q_all.txt",0)' </dev/null 2>&1 | tail -1
echo QALL_DONE
