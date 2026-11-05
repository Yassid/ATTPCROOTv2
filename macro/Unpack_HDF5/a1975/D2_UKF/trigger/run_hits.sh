#!/bin/bash
source /mnt/f/a1975/sim_rebuild/r280/env.sh; cd /mnt/f/a1975/sim_rebuild/trig; mkdir -p qhits
root -l -b -q dump_sim_tracks.C </dev/null >/dev/null 2>&1
awk '$1==0{print $2}' sim_tracks.txt | sort -n -u > qhits/ev_s0.txt; awk '$1==1{print $2}' sim_tracks.txt | sort -n -u > qhits/ev_s1.txt
awk '$3==1 && $14==1{print $1}' tracks.txt | sort -n -u > qhits/runs.txt
for r in $(cat qhits/runs.txt); do awk -v r=$r '$1==r && $3==1 && $14==1{print $2}' tracks.txt | sort -n -u > qhits/ev_$r.txt; done
{ echo "/mnt/f/a1975/sim_rebuild/r280/p_gs_reco.root qhits/ev_s0.txt qhits/q_s0.txt 0"; echo "/mnt/f/a1975/sim_rebuild/r280/p_ex39_reco.root qhits/ev_s1.txt qhits/q_s1.txt 1"; for r in $(cat qhits/runs.txt); do printf "/mnt/f/a1975/reco_d2_dv1104/run_%04d_multifit_reco.root qhits/ev_%d.txt qhits/q_%d.txt %d\n" $r $r $r $r; done; } > qhits/jobs.txt
cat qhits/jobs.txt | xargs -P 6 -L 1 bash -c 'root -l -b -q "mesh_q_hits.C(\"$0\",\"$1\",\"$2\",$3)" </dev/null 2>&1 | grep -E "^done" '
echo HITS_ALLDONE
