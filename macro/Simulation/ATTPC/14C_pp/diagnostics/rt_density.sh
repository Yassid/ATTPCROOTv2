#!/usr/bin/env bash
# Gas-density test on the corrected (negB) sample, UKF only.
#   3.553e-5 = the 273 K value currently used everywhere (media.geo H_300torr)
#   3.308e-5 = the 293 K "RT" value; the matFX study measured the former as ~7 % too dense
# UKF is the right fitter for this: its GetKinematics() is back-extrapolated to the vertex using
# dEdx*pathLength from the SAME density, so the test exercises transport and vertex correction
# together. (GENFIT's GetKinematics() is the uncorrected first-point value -- separate issue.)
set -eo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); SIM=$(cd "$HERE/.." && pwd)
REPO=$(cd "$SIM/../../../.." && pwd); UKF="$REPO/macro/Unpack_HDF5/a1954/UKF"
IN="$HERE/negB/"; NAME=negBg
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
cd "$SIM"
ln -sf "${IN}${NAME}_reco.root" "${IN}${NAME}RT_reco.root"
echo "[$(date +%H:%M:%S)] UKF, gasDensity = 3.308e-5 (RT)"
root -b -q -l "$UKF/pipeline/fitUKF_C14.C(\"${NAME}RT\",-1,\"proton\",-1,2.85,3.308e-5,\"\",\"$IN\",0.5,0.1,1,10,\"$IN\")" 2>&1 | tail -2
root -b -q -l "$UKF/pp/ex_C14.C(\"${NAME}RT\",\"$IN\",161.0,1e9,\"_${NAME}RT\",1.007825,14.003242,\"\",\"ukf\")" 2>&1 | grep "good track"
echo RT_COMPLETED
