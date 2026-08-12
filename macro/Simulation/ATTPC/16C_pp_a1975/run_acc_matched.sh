#!/usr/bin/env bash
# Reco + fit + acceptance with the chain MATCHED to the data:
#   AtPSAMax(20), no cleaning, AtTrackFinderTC(15.0, 7.5), proton PID gate in the fitter.
# Space charge is deliberately not reproduced.
#
# The acceptance uses the beam energy the sample was GENERATED at (192), not the one the data
# analysis assumes (188). Mixing them puts truth events in the wrong theta_cm bin.
set -eo pipefail
SIM=${1:-/mnt/f/a1975_C16_pp_acc/el_s1002_sim.root}
OUT=${2:-/mnt/f/a1975_C16_pp_matched}
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
mkdir -p "$OUT"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"

echo "[$(date +%H:%M:%S)] reco  (AtPSAMax, no cleaning, TC)"
cd "$REPO/macro/Simulation/ATTPC/16C_pp_a1975"
root -b -q -l "run_reco_C16_TC.C(\"$SIM\",\"$OUT/m_reco.root\",\"ATTPC.a1975_C16_sim.par\",20,0,20.0)" > "$OUT/reco.log" 2>&1
grep -q "sim reco done" "$OUT/reco.log" || { echo RECO_FAILED; exit 1; }

echo "[$(date +%H:%M:%S)] fit   (proton hypothesis + proton_band gate)"
cd "$REPO/macro/Unpack_HDF5/a1975/UKF"
root -b -q -l "pipeline/fitGenfitter_a1975.C(\"m\",-1,\"$OUT/\",\"\",\"$OUT/\",-2.85,2,5,\"proton_H2_catima.txt\",\"pid/proton_band.json\",4.0,kFALSE,2212,1.00727647,1,kFALSE,\"ATTPC_H300torr_RT\")" > "$OUT/fit.log" 2>&1
[ -s "$OUT/m_genfitter.root" ] || { echo FIT_FAILED; exit 1; }

echo "[$(date +%H:%M:%S)] acceptance"
cd "$REPO/macro/Simulation/ATTPC/16C_pp_a1975"
root -b -q -l "acceptance_C16pp.C(\"$SIM\",\"$OUT/m_genfitter.root\",\"matched\",0.0,192.0,5.0,36,180.0,10.0,0.5,2.0,kFALSE)" > "$OUT/acc.log" 2>&1
grep "overall acceptance" "$OUT/acc.log"
rm -f "$OUT/m_reco.root"
echo "[$(date +%H:%M:%S)] MATCHEDDONE"
