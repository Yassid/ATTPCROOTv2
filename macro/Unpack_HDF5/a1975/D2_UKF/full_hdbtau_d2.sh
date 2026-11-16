#!/usr/bin/env bash
# 2026-11-16: 16C(d,p) on the D2 runs re-unpacked with HDBSCAN + float tau (reco_hdbtau_d2_local.sh), then the ADOPTED
# (d,p) genfit (find805: FIND ladder 80 -> 5 % in 5 % steps, MINCL 8, C2MAX 0.1, peak-amplitude proton gate
# pid/proton_d2_dv1104.json, CATIMA, r280 fit geometry as gf_dp_find805_r280) into gf_dp_hdbtau/, then the (d,p) cache dp_kin_hdbtau.root.
# Resume = rerun (every stage skips what is done). Lock as full_hdbtau.sh.
D=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975
mkdir -p /home/yassid/d2stage /mnt/f/a1975/logs_d2_hdbtau
exec 9>/home/yassid/d2stage/full_hdbtau_d2.lock
flock -n 9 || { echo "full_hdbtau_d2.sh: another instance (or its orphaned workers) holds the lock -- abort"; exit 1; }
cd $D/D2_UKF   # 2026-11-16 copy of full_tctau_d2.sh: TC replaced by HDBSCAN (TC fragments backward proton spirals)
NPAR=6 MAXSTAGED=6 ./reco_hdbtau_d2_local.sh &
sleep 5; GPAT='reco_hdbtau_d2_local[.]sh' GLOG=/mnt/f/a1975/logs_d2_hdbtau/guard_cfree.log ./guard_cfree_tc.sh &
wait
echo "=== reco stage over $(date '+%F %T'): $(ls /mnt/f/a1975/reco_d2_hdbtau/*.taudone 2>/dev/null | wc -l) / 47 ==="
REC=/mnt/f/a1975/reco_d2_hdbtau/ GF=/mnt/f/a1975/gf_dp_hdbtau/ LOG=/mnt/f/a1975/logs_dp_hdbtau/ GEO=ATTPC_D300torr_r280_geomanager.root \
  MAXPCT=80 MINPCT=5 STEP=5 ./prod_dp_findscan.sh 6
mkdir -p /home/yassid/d2stage/cache_work_hdb
FIT=/mnt/f/a1975/gf_dp_hdbtau/ OUT=/mnt/f/a1975/caches/dp_kin_hdbtau.root ./cache_dp_dv1104.sh 3 /home/yassid/d2stage/cache_work_hdb
echo "=== full_hdbtau_d2 finished $(date '+%F %T') ==="
