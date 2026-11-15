#!/usr/bin/env bash
# 2026-11-14: FULL TC + float-tau + integral-gate production, runs 0106-0189 (the 3 test runs are skipped as done):
# reco (NPAR 6, MAXSTAGED 4, C: guard) -> (p,d) genfit with the integral gate -> event-identity caches.
D=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975
# 2026-11-14: refuse a second instance. Every child (xargs workers included) inherits fd 9, so the lock
# stays held while ANY process of a previous launch is alive -- an orphaned worker set from a killed
# launch reco'd every run twice into the same output file.
exec 9>/home/yassid/h2stage/full_tctau.lock
flock -n 9 || { echo "full_tctau.sh: another instance (or its orphaned workers) holds the lock -- abort"; exit 1; }
mkdir -p /mnt/f/a1975/logs_h2_tctau
cd $D/D2_UKF
NPAR=6 MAXSTAGED=6 ./reco_tctau_h2_local.sh &
sleep 5; ./guard_cfree_tc.sh &
wait
cd $D/UKF/pd
NPAR=4 ./gf_pd_h2tcint.sh
cd /home/yassid/a1975_analysis/pd/production
./cache_pd_tcint.sh
echo "=== full_tctau finished $(date '+%F %T') ==="
