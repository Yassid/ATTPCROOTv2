#!/usr/bin/env bash
# GENFIT on the corrected (negB) sample WITH the CATIMA vertex-gap back-extrapolation, the
# a1975 D2 pattern: matEffects stays OFF (keeps genfit's ~85 % good-fit rate) while CATIMA
# supplies the dE/dx over the geometric vertex gap. Result lands in GetKinematicsXtr(), so the
# Ex step must be told to read that slot (useXtr = kTRUE) -- reading GetKinematics() would give
# the uncorrected first-cluster energy, which is what every previous 14C genfit number used.
# H2 target -> matA = 1. Density at RT = 3.308e-5 g/cm3.
set -eo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); SIM=$(cd "$HERE/.." && pwd)
REPO=$(cd "$SIM/../../../.." && pwd); UKF="$REPO/macro/Unpack_HDF5/a1954/UKF"
IN="$HERE/negB/"; RHO=${1:-3.308e-5}; TAG=${2:-negBgXtr}
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
cd "$SIM"
ln -sf "${IN}negBg_reco.root" "${IN}${TAG}_reco.root"
echo "[$(date +%H:%M:%S)] GENFIT + backExtrap + CATIMA (rho=$RHO, A=1), matEffects OFF"
root -b -q -l "$UKF/pipeline/fitGenfit_C14.C(\"$TAG\",-1,\"$IN\",\"\",\"$IN\",-2.85,2,5,\"\",4.0,10.0,170.0,kFALSE,kFALSE,\"proton\",\"ATTPC_H300torr\",kTRUE,kTRUE,$RHO,1)" 2>&1 | grep -E "Back-extrap|CATIMA|Done|Real" | head -4
# both slots, so the size of the correction is measurable rather than assumed
root -b -q -l "$UKF/pp/ex_C14.C(\"$TAG\",\"$IN\",161.0,1e9,\"_${TAG}_raw\",1.007825,14.003242,\"\",\"genfit\",kFALSE)" 2>&1 | grep "good track"
root -b -q -l "$UKF/pp/ex_C14.C(\"$TAG\",\"$IN\",161.0,1e9,\"_${TAG}_xtr\",1.007825,14.003242,\"\",\"genfit\",kTRUE)"  2>&1 | grep "good track"
echo CATIMA_COMPLETED
