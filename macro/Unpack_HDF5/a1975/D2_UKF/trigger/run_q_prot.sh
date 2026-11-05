#!/bin/bash
source /mnt/f/a1975/sim_rebuild/r280/env.sh; cd /mnt/f/a1975/sim_rebuild/trig
root -l -b -q 'mesh_q_prot.C("/mnt/f/a1975/sim_rebuild/r280/p_gs_reco.root","qhits/ev_all.txt","qhits/qp_all.txt",0)' </dev/null 2>&1 | tail -1 &
root -l -b -q 'mesh_q_prot.C("/mnt/f/a1975/sim_rebuild/r280/p_ex39_reco.root","qhits/ev_s1.txt","qhits/qp_s1.txt",1)' </dev/null 2>&1 | tail -1 &
wait; echo QPROT_DONE
