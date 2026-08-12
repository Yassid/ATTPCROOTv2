#!/usr/bin/env bash
# Full a1975 (p,p) chain with the DATA's track finder (AtTrackFinderTC) and the DATA's PID gate.
#
# Both were mismatched before: the simulation used AtTrackFinderHDBSCAN (mcs 20), which discards
# the short tracks the forward-angle recoil protons make, and the fit ran with no PID gate while
# the data production gates on pid/proton_band.json inside the fitter. Neither the acceptance nor
# any comparison of track populations means anything until both match.
set -eo pipefail
SIM=${1:-/mnt/f/a1975_C16_pp_acc/el_s1002_sim.root}
OUT=${2:-/mnt/f/a1975_C16_pp_tc}
GATE=${3:-pid/proton_band.json}
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
mkdir -p "$OUT"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"

echo "[$(date +%H:%M:%S)] reco with AtTrackFinderTC"
cd "$REPO/macro/Simulation/ATTPC/16C_pp_a1975"
root -b -q -l "run_reco_C16_TC.C(\"$SIM\",\"$OUT/tc_reco.root\",\"ATTPC.a1975_C16_sim.par\",20,0,20.0)" > "$OUT/reco.log" 2>&1
grep -q "sim reco done" "$OUT/reco.log" || { echo "RECO_FAILED"; exit 1; }

cd "$REPO/macro/Unpack_HDF5/a1975/UKF"
# ungated AND gated, so the gate's effect is separable from the track finder's
for G in "" "$GATE"; do
  TAG=$([ -z "$G" ] && echo "nogate" || echo "gated")
  echo "[$(date +%H:%M:%S)] fit ($TAG)"
  root -b -q -l "pipeline/fitGenfitter_a1975.C(\"tc\",-1,\"$OUT/\",\"_$TAG\",\"$OUT/\",-2.85,2,5,\"proton_H2_catima.txt\",\"$G\",4.0,kFALSE,2212,1.00727647,1,kFALSE,\"ATTPC_H300torr_RT\")" > "$OUT/fit_$TAG.log" 2>&1
  [ -s "$OUT/tc_genfitter_$TAG.root" ] || echo "  $TAG FIT_FAILED"
done
echo "[$(date +%H:%M:%S)] TCDONE"
