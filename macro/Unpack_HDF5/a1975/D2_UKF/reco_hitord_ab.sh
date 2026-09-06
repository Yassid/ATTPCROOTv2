#!/usr/bin/env bash
# Re-reco of the a1975 D2 runs with AtPRA::OrderHitsAlongTrack ACTIVE in the HDBSCAN path.
#
# EVERY argument is dv_dt.sh's dv1104 production line verbatim -- thr 40, hdbscan, mcs 0, ms 3,
# par ATTPC.a1975_deuterium_dv1104.par -- so the ONLY difference from /mnt/f/a1975/reco_d2_dv1104/
# is the hit ordering before SetTrackInitialParameters. Output goes to a SEPARATE directory;
# production is never overwritten, so the old file is the control for the A/B.
#
#   ./reco_hitord_ab.sh [runs...]
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
set -uo pipefail

RUNS="${*:-0020 0016}"
H5=/mnt/f/a1975/h5/
REC=/mnt/f/a1975/reco_d2_dv1104_hitord/
LOG=/mnt/f/a1975/logs_hitord/
PAR=ATTPC.a1975_deuterium_dv1104.par
mkdir -p "$REC" "$LOG"

one() {
  r="run_$1"
  rc="${REC}${r}_multifit_reco.root"
  # resume on a marker that is only written after root exits AND the product reads back non-trivial
  [ -s "${rc}.hitorddone" ] && { echo "[have] $r"; return 0; }
  echo "[start] $r $(date '+%H:%M:%S')"
  root -l -b -q "unpackReco_multifit.C(\"$r\",0,false,\"$REC\",\"$H5\",false,true,\"multifit\",0,40,\"hdbscan\",0,3,0,0.1,\"$PAR\")" \
      > "${LOG}reco_${r}.log" 2>&1
  if grep -q 'segmentation violation' "${LOG}reco_${r}.log" || [ ! -s "$rc" ] \
     || [ "$(stat -c%s "$rc" 2>/dev/null || echo 0)" -lt 20000000 ]; then
    echo "[FAIL] $r  (see ${LOG}reco_${r}.log)"; rm -f "$rc"; return 1
  fi
  touch "${rc}.hitorddone"
  echo "[done ] $r $(date '+%H:%M:%S')  $(du -h "$rc" | cut -f1)"
}
export -f one; export REC LOG PAR H5

printf '%s\n' $RUNS | xargs -P 2 -I{} bash -c 'one "$@"' _ {}
echo "=== reco_hitord_ab finished $(date '+%H:%M:%S') ==="
