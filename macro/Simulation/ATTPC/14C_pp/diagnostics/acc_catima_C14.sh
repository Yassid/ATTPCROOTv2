#!/usr/bin/env bash
# 14C acceptance regenerated against the 2026-08-25 CATIMA production. One (level, seed) per call.
#
# WHAT CHANGES vs acc_batch_genfit.sh, and why the truth must be REGENERATED rather than reused:
#   * transport gas ATTPC_H300torr_RT (3.308e-5), not ATTPC_H300torr (3.553e-5, a 0 C value).
#     The old truth was generated with 7.4 % too much material, so it is not reusable -- the
#     energy loss is in the truth, not just in the fit.
#   * genfit with matEffects = kTRUE + native CATIMA dE/dx, matFallback = kFALSE, geometry
#     ATTPC_H300torr_RT, and NO manual gap eloss (extrapolateToLine already integrates the gap
#     once material effects are on, so the pair double-counts it).
# The acceptance itself is still evaluated at chi2/ndf < 5 with useXtr = kTRUE. The chi2 choice
# barely matters: the chi2-cut and no-cut acceptances agree to 0.2-0.3 % except at 22 deg.
#
#   ./acc_catima_C14.sh <tag> <resEx> <seed> [nEvents]
set -eo pipefail
TAG=$1; EX=$2; SEED=$3; NEV=${4:-8000}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); SIM=$(cd "$HERE/.." && pwd)
REPO=$(cd "$SIM/../../../.." && pwd); UKF="$REPO/macro/Unpack_HDF5/a1954/UKF"
OUT=${ACC_OUT:-/mnt/f/a1954_C14_acc_catima}; mkdir -p "$OUT"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
cd "$SIM"
J="${TAG}_s${SEED}"
[ -f "$OUT/$J.marker" ] && { echo "$J already COMPLETED"; exit 0; }
if [ ! -s "$OUT/${J}_sim.root" ]; then
  root -b -q -l "C14_pp_sim.C($NEV,2.0,178.0,\"TGeant4\",-28.5,\"$OUT/${J}_sim.root\",$EX,$SEED,\"ATTPC_H300torr_RT.root\")" > "$OUT/${J}_gen.log" 2>&1
  grep -q "ATTPC_H300torr_RT" "$OUT/${J}_gen.log" || { echo "$J WRONG_GAS"; exit 1; }
  grep -q "RNG seed requested: $SEED" "$OUT/${J}_gen.log" || { echo "$J SEED_NOT_APPLIED"; exit 1; }
fi
[ -s "$OUT/${J}_reco.root" ] || root -b -q -l "run_reco_C14.C(\"$OUT/${J}_sim.root\",\"$OUT/${J}_reco.root\",\"ATTPC.a1954_C14_sim.par\",20,20,8,0,30.0,\"mover\")" > "$OUT/${J}_reco.log" 2>&1
[ -s "$OUT/${J}_reco.root" ] || { echo "$J RECO_FAILED"; exit 1; }
root -b -q -l "$UKF/pipeline/fitGenfit_C14.C(\"$J\",-1,\"$OUT/\",\"\",\"$OUT/\",-2.85,2,5,\"\",4.0,10.0,170.0,kTRUE,kFALSE,\"proton\",\"ATTPC_H300torr_RT\",kFALSE,kTRUE,0.0,1,kTRUE,kFALSE,kFALSE,kFALSE)" > "$OUT/${J}_fit.log" 2>&1
[ -s "$OUT/${J}_genfit.root" ] || { echo "$J FIT_FAILED"; exit 1; }
grep -q "dE/dx from CATIMA" "$OUT/${J}_fit.log" || { echo "$J CATIMA_NOT_ENABLED"; exit 1; }
root -b -q -l "acceptance_C14.C(\"$OUT/${J}_sim.root\",\"$OUT/${J}_genfit.root\",\"$J\",$EX,159.75,5.0,36,180.0,10.0,0.5,2.0,kTRUE)" > "$OUT/${J}_acc.log" 2>&1
grep -q "overall acceptance" "$OUT/${J}_acc.log" || { echo "$J ACC_FAILED"; exit 1; }
mv -f "$HERE/acceptance_${J}.root" "$OUT/" 2>/dev/null || true
rm -f "$HERE/acceptance_${J}.png" "$OUT/${J}_reco.root"
echo COMPLETED > "$OUT/$J.marker"
echo "$J done: $(grep 'overall acceptance' "$OUT/${J}_acc.log")"
