#!/usr/bin/env bash
# 16C(d,t)15C at DriftVelocity = 1.35 cm/us.
#
# Reco: same production config as reco_d2_hdb1136 (thr 40, hdbscan, mcs 0, ms 3) so the ONLY
# difference from that set is the drift velocity -- the par is ATTPC.a1975_deuterium_dv${DVTAG}.par
# and it is passed to BOTH the reco and the fit (the fit macro used to hardcode the 1.136 par).
#
# Fit: genfit with matEffects OFF, back-extrapolation to the beam axis ON, and the hand-applied
# CATIMA correction over the vertex gap (D2 at 300 torr, rho 6.61e-5, A=2). With material effects
# off the extrapolation is geometric and leaves |p| untouched, so CATIMA is the only thing that
# makes the vertex energy differ from the fitted one -- which is the comparison being set up:
#   ke    (AtFittedTrack Xtr slot) = back-extrapolated + CATIMA
#   kefit (AtFittedTrack raw slot) = what the fit returned at the first measurement point
#
#   ./dv${DVTAG}_dt.sh [nparallel] [runs...]
DVTAG="${DVTAG:-135}"
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
# config.sh trips over unset vars, so it must be sourced BEFORE set -u (it silently kills the
# whole script otherwise). It also omits ROOT_INCLUDE_PATH; without build/include on it, cling
# crashes on the hidden-friend swap when a macro loads the libraries.
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
set -u

NPAR="${1:-4}"; shift || true
RUNS="${*:-0016 0017 0018 0019 0020 0021 0022 0023 0026 0027 0031 0032 0034 0036 0037 0038 0039 0040 0041 0042 0043 0044 0046 0048 0057 0058 0076 0077 0078 0079 0080 0082 0083 0084 0085 0086 0087 0088 0089 0091 0092 0095 0096 0097 0098 0102 0103}"

H5=/mnt/f/a1975/h5/
REC=/mnt/f/a1975/reco_d2_dv${DVTAG}/
FIT=/mnt/f/a1975/gf_dt_dv${DVTAG}/
LOG=/mnt/f/a1975/logs_dv${DVTAG}/
PAR=ATTPC.a1975_deuterium_dv${DVTAG}.par
mkdir -p "$REC" "$FIT" "$LOG"

one() {
  n="$1"; r="run_${n}"
  REC=/mnt/f/a1975/reco_d2_dv${DVTAG}/; FIT=/mnt/f/a1975/gf_dt_dv${DVTAG}/; LOG=/mnt/f/a1975/logs_dv${DVTAG}/
  PAR=ATTPC.a1975_deuterium_dv${DVTAG}.par
  rc="${REC}${r}_multifit_reco.root"; fo="${FIT}${r}_multifit_genfitter_t.root"

  # A marker means COMPLETED, never "a file exists": the macro's own Done line, a size floor,
  # and no segfault in the log. A failed stage deletes its own partial output.
  if [ ! -s "${rc}.dv${DVTAG}done" ]; then
    root -l -b -q "unpackReco_multifit.C(\"$r\",0,false,\"$REC\",\"$H5\",false,true,\"multifit\",0,40,\"hdbscan\",0,3,0,0.1,\"$PAR\")" \
        > "${LOG}reco_${r}.log" 2>&1
    if grep -q 'segmentation violation' "${LOG}reco_${r}.log" || [ ! -s "$rc" ] \
       || [ "$(stat -c%s "$rc" 2>/dev/null || echo 0)" -lt 20000000 ]; then
      echo "[FAIL reco] $r"; rm -f "$rc"; return 1
    fi
    touch "${rc}.dv${DVTAG}done"
  fi

  if [ ! -s "${fo}.dv${DVTAG}done" ]; then
    root -l -b -q "fitGenfitter_a1975_deuterium.C(\"${r}_multifit\",-1,\"$REC\",\"\",\"$FIT\",-2.85,2,5,\"pid/triton_d2.json\",4.0,10.0,170.0,kFALSE,kTRUE,1000010030,3.01550072,1,\"t\",\"_reco\",\"ATTPC_D300torr_v2_geomanager.root\",kTRUE,6.61e-5,2,\"$PAR\")" \
        > "${LOG}fit_${r}.log" 2>&1
    if grep -q 'segmentation violation' "${LOG}fit_${r}.log" || [ ! -s "$fo" ]; then
      echo "[FAIL fit] $r"; rm -f "$fo"; return 1
    fi
    touch "${fo}.dv${DVTAG}done"
  fi
  echo "[done] $r  $(date '+%H:%M:%S')"
}
export -f one; export H5 DVTAG
printf '%s\n' $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "=== dv${DVTAG}: $(ls ${FIT}*_genfitter_t.root 2>/dev/null | wc -l) fits in $FIT ==="
