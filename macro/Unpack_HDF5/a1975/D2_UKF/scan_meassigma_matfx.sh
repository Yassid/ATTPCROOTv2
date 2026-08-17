#!/usr/bin/env bash
# Does the matFX fit collapse depend on the assumed per-cluster position error?
#
# AtGenfitter builds the measurement covariance as
#   varT = measSigma^2 + fDiffTransMM^2 * L_drift_cm
#   varZ = fZLongFactor*measSigma^2 + fDiffLongMM^2 * L_drift_cm
# and the diffusion coefficients DEFAULT TO ZERO, so in production the covariance is flat
# (s2, s2, 2*s2) for every cluster no matter how far its electrons drifted. If the collapse is
# the filter rejecting measurements because it was handed over-tight errors, the collapse rate
# must move as measSigma is varied -- and must move DIFFERENTLY for a run that collapses (0016,
# 85%) and one that does not (0031, 3%).
#
# A flat response in both runs kills that hypothesis too.
#
# Everything except measSigma is the production configuration: matFX on, CATIMA MSC +
# straggling, the corrected dE/dx table, manualElossDensity 0 (the gap loss is already in the
# RK back-extrapolation), matFallback off so a throwing track is dropped rather than silently
# refitted without material.
#
#   ./scan_meassigma_matfx.sh [nEvents] [nparallel]
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
# config.sh reads unset vars: source it BEFORE set -u or it kills the script silently
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export GENFIT=/home/yassid/fair_install/GenFit
export LD_LIBRARY_PATH=$GENFIT/lib:${LD_LIBRARY_PATH:-}
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
set -uo pipefail

NEV="${1:-5000}"
NPAR="${2:-6}"
SIGMAS="1.5 2.5 4.0 6.0 9.0 14.0"
RUNS="0016 0031"
REC=/mnt/f/a1975/reco_d2_dv1104/
PAR=ATTPC.a1975_deuterium_dv1104.par
GATE=pid/triton_d2_dv1104.json
TAB=triton_D2_300torr.txt
BASE=/mnt/f/a1975/sigscan
mkdir -p "$BASE"

one() {
  r="$1"; s="$2"
  tag="run_${r}_s${s}"
  out="${BASE}/${tag}/"
  mkdir -p "$out"
  f="${out}run_${r}_multifit_genfitter_t.root"
  [ -s "$f" ] && { echo "[have] $tag"; return 0; }
  root -l -b -q "fitGenfitter_a1975_deuterium.C(\"run_${r}_multifit\",${NEV},\"$REC\",\"\",\"$out\",\
-2.85,2,5,\"$GATE\",${s},10.0,170.0,kTRUE,kTRUE,1000010030,3.01550072,1,\"t\",\"_reco\",\
\"ATTPC_D300torr_v2_geomanager.root\",kTRUE,0,2,\"$PAR\",kFALSE,kFALSE,kFALSE,\"$TAB\",kTRUE,kTRUE)" \
    > "${BASE}/${tag}.log" 2>&1
  [ -s "$f" ] && echo "[ok] $tag" || echo "[FAIL] $tag"
}
export -f one
export BASE NEV REC PAR GATE TAB

echo "=== measSigma scan: runs $RUNS x sigma {$SIGMAS}, $NEV events each ==="
for r in $RUNS; do for s in $SIGMAS; do echo "$r $s"; done; done \
  | xargs -P "$NPAR" -n 2 bash -c 'one "$0" "$1"'
echo "=== fits done; tallying ==="
root -l -b -q "tally_sigscan.C(\"$BASE\",\"$SIGMAS\",\"$RUNS\")"
