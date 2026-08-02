#!/usr/bin/env bash
# Per-run UKF fit worker (D2 (d,t), triton hyp) -> reco_d2/<run>_genfitter_t_UKF.root
# Input is the multifit point cloud (_multifit_reco), gated on its .mfdone marker.
set -e
RUN="$1"
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
IODIR=/mnt/f/a1975/reco_d2/
OUT="${IODIR}${RUN}_genfitter_t_UKF.root"
if [ ! -f "${IODIR}${RUN}_multifit_reco.root.mfdone" ]; then echo "[wait] $RUN multifit reco not done"; exit 0; fi
if [ -f "${OUT}.done" ]; then echo "[skip] $RUN ukf-t (done)"; exit 0; fi
# A .root with no .done is a run that was interrupted mid-write (e.g. the 2026-08-02
# crash): ROOT can only open it in recovery mode with a partial tree. Discard and redo.
if [ -f "$OUT" ]; then echo "[redo] $RUN ukf-t (truncated leftover removed)"; rm -f "$OUT"; fi
source "$REPO/build/config.sh" >/dev/null 2>&1
cd "$REPO/macro/Unpack_HDF5/a1975/D2_UKF"
echo "[start] ukf-t $RUN $(date '+%H:%M:%S')"
root -l -b -q "fitUKF_a1975_deuterium.C(\"$RUN\", -1, \"triton\", -1, 2.85, 9.0e-5, \"$IODIR\", \"$IODIR\", 2.0, 0.3, 1, 4, kTRUE, \"_multifit_reco\", \"_genfitter_t_UKF\")" \
   > "ukflog_t_${RUN}.txt" 2>&1
touch "${OUT}.done"
echo "[done ] ukf-t $RUN $(date '+%H:%M:%S') -> $(du -h "$OUT" 2>/dev/null|cut -f1)"
