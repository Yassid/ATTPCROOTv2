#!/usr/bin/env bash
# GENFIT acceptance, one job = one level + one RNG seed.
#
# Same MC truth as the UKF acceptance (the *_sim.root files are reused), so the two acceptances
# differ ONLY by the fitter and can be compared directly. The reco has to be remade because the
# UKF pass deleted it to save space; it is deterministic given the same sim + par, so the point
# cloud is identical to the one the UKF fits saw.
#
# GENFIT runs with matEffects OFF plus SetBackExtrapToAxis + CATIMA manual eloss over the vertex
# gap (rho = 3.308e-5 g/cm3, matA = 1 for the H2 target). That is the a1975 D2 pattern: it keeps
# genfit's ~85 % good-fit rate instead of the ~60 % that material effects cost, while still
# recovering the energy lost before the first cluster. The corrected energy lands in
# GetKinematicsXtr(), so the acceptance is evaluated with useXtr = kTRUE -- genfit's plain
# GetKinematics() is the RAW first-cluster value and would be the wrong quantity.
#
#   ./acc_batch_genfit.sh <tag> <resEx> <seed> [nEvents]
set -eo pipefail
TAG=$1; EX=$2; SEED=$3; NEV=${4:-8000}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); SIM=$(cd "$HERE/.." && pwd)
REPO=$(cd "$SIM/../../../.." && pwd); UKF="$REPO/macro/Unpack_HDF5/a1954/UKF"
SRC=/mnt/f/a1954_C14_acc            # where the shared MC truth lives
OUT=${ACC_OUT:-/mnt/f/a1954_C14_acc_gf}
mkdir -p "$OUT"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
cd "$SIM"
J="${TAG}_s${SEED}"
[ -f "$OUT/$J.marker" ] && { echo "$J already COMPLETED"; exit 0; }

[ -s "$SRC/${J}_sim.root" ] || { echo "$J NO_TRUTH ($SRC/${J}_sim.root)"; exit 1; }
ln -sf "$SRC/${J}_sim.root" "$OUT/${J}_sim.root"

if [ ! -s "$OUT/${J}_reco.root" ]; then
  root -b -q -l "run_reco_C14.C(\"$OUT/${J}_sim.root\",\"$OUT/${J}_reco.root\",\"ATTPC.a1954_C14_sim.par\",20,20,8,0,30.0,\"mover\")" > "$OUT/${J}_reco.log" 2>&1
  [ -s "$OUT/${J}_reco.root" ] || { echo "$J RECO_FAILED"; exit 1; }
fi

# fitGenfit_C14.C(name,nEv,ioDir,outSuffix,outDir,bField,minIter,maxIter,pidGate,measSigma,
#                 thMin,thMax,matEffects,backwardSeedFix,particle,geoName,matFallback,
#                 backExtrap,manualElossDensity,matA)
root -b -q -l "$UKF/pipeline/fitGenfit_C14.C(\"$J\",-1,\"$OUT/\",\"\",\"$OUT/\",-2.85,2,5,\"\",4.0,10.0,170.0,kFALSE,kFALSE,\"proton\",\"ATTPC_H300torr\",kTRUE,kTRUE,3.308e-5,1)" > "$OUT/${J}_fit.log" 2>&1
[ -s "$OUT/${J}_genfit.root" ] || { echo "$J FIT_FAILED"; exit 1; }
grep -q "CATIMA" "$OUT/${J}_fit.log" || { echo "$J CATIMA_NOT_ENABLED"; exit 1; }

root -b -q -l "acceptance_C14.C(\"$OUT/${J}_sim.root\",\"$OUT/${J}_genfit.root\",\"$J\",$EX,161.0,5.0,36,180.0,10.0,0.5,2.0,kTRUE)" > "$OUT/${J}_acc.log" 2>&1
grep -q "overall acceptance" "$OUT/${J}_acc.log" || { echo "$J ACC_FAILED"; exit 1; }
mv -f "$HERE/acceptance_${J}.root" "$OUT/" 2>/dev/null || true
rm -f "$HERE/acceptance_${J}.png" "$OUT/${J}_reco.root"
echo COMPLETED > "$OUT/$J.marker"
echo "$J done: $(grep 'overall acceptance' "$OUT/${J}_acc.log")"
