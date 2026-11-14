#!/usr/bin/env bash
# 2026-11-14: TC + float-tau + integral-gate TEST on 3 runs (0107 0109 0113): reco, then (p,d) genfit.
D=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975
mkdir -p /mnt/f/a1975/logs_h2_tctau
RUNS="0107 0109 0113"
cd $D/D2_UKF
NPAR=3 MAXSTAGED=3 ./reco_tctau_h2_local.sh $RUNS &
sleep 5; ./guard_cfree_tc.sh &
wait
cd $D/UKF/pd
NPAR=3 ./gf_pd_h2tcint.sh $RUNS
echo "=== test3_tctau finished $(date '+%F %T') ==="
