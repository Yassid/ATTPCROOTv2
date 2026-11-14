#!/usr/bin/env bash
# d/t test of the trace-integral PID (C15d option) on more a1975 D2 runs.
#
# Re-reco with the pulse peaking time FLOATED (fitChi2 = -1: float tau, no chi2 gate), every other
# argument identical to the run_0016 test in /mnt/f/a1975/reco_d2_tau/ (and to the dv1104
# production apart from fitChi2). The production (fitChi2 = 0) has a fixed tau, so its integral is
# 2.9618 x amplitude and cannot be used. Then pid_width_dp.C writes the per-track PID tree with
# sqa (peak amplitude) and sqw (width-weighted = integral) for the same tracks.
#
#   ./reco_tau_d2.sh [runs...]      output /mnt/f/a1975/reco_d2_tau/, logs /mnt/f/a1975/logs_tau/
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
set -uo pipefail

RUNS="${*:-0017 0018 0019 0020 0021 0022 0024 0025 0026 0027 0031 0032}"
H5=/mnt/f/a1975/h5/
REC=/mnt/f/a1975/reco_d2_tau/
LOG=/mnt/f/a1975/logs_tau/
PAR=ATTPC.a1975_deuterium_dv1104.par
NPAR=${NPAR:-6}
mkdir -p "$REC" "$LOG"

one() {
  r="run_$1"
  rc="${REC}${r}_multifit_reco.root"
  pw="${REC}pid_width_$1.root"
  [ -s "${pw}.done" ] && { echo "[have] $r"; return 0; }
  if [ ! -s "${rc}.taudone" ]; then
    echo "[reco ] $r $(date '+%H:%M:%S')"
    root -l -b -q "unpackReco_multifit.C(\"$r\",0,false,\"$REC\",\"$H5\",false,true,\"multifit\",0,40,\"hdbscan\",0,3,-1,0.1,\"$PAR\")" \
        > "${LOG}reco_${r}.log" 2>&1
    if grep -q 'segmentation violation' "${LOG}reco_${r}.log" || [ ! -s "$rc" ] \
       || [ "$(stat -c%s "$rc" 2>/dev/null || echo 0)" -lt 20000000 ]; then
      echo "[FAIL reco] $r"; rm -f "$rc"; return 1
    fi
    touch "${rc}.taudone"
  fi
  echo "[pid  ] $r $(date '+%H:%M:%S')"
  root -l -b -q "pid_width_dp.C(\"$rc\",\"$pw\")" > "${LOG}pid_${r}.log" 2>&1
  if [ ! -s "$pw" ]; then echo "[FAIL pid] $r"; return 1; fi
  touch "${pw}.done"
  echo "[done ] $r $(date '+%H:%M:%S')  $(du -h "$rc" | cut -f1)"
}
export -f one; export REC LOG PAR H5

printf '%s\n' $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "=== reco_tau_d2 finished $(date '+%H:%M:%S') ==="
