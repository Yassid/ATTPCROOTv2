#!/usr/bin/env bash
# Per-run UKF fit worker (D2 (d,p), proton hyp) -> reco_d2/<run>_genfitter_p_UKF.root
set -e
RUN="$1"
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
IODIR=/mnt/f/a1975/reco_d2/
OUT="${IODIR}${RUN}_genfitter_p_UKF.root"
if [ ! -f "${IODIR}${RUN}_reco.root.done" ]; then echo "[wait] $RUN reco not done"; exit 0; fi
if [ -f "${OUT}.done" ]; then echo "[skip] $RUN ukf (done)"; exit 0; fi
source "$REPO/build/config.sh" >/dev/null 2>&1
cd "$REPO/macro/Unpack_HDF5/a1975/D2_UKF"
echo "[start] ukf $RUN $(date '+%H:%M:%S')"
root -l -b -q "fitUKF_a1975_deuterium.C(\"$RUN\", -1, \"proton\", -1, 2.85, 9.0e-5, \"$IODIR\")" \
   > "ukflog_${RUN}.txt" 2>&1
touch "${OUT}.done"
echo "[done ] ukf $RUN $(date '+%H:%M:%S') -> $(du -h "$OUT" 2>/dev/null|cut -f1)"
