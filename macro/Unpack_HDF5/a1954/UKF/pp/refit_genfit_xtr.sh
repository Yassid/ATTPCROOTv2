#!/usr/bin/env bash
# Refit the a1954 14C data with GENFIT in EXACTLY the configuration the GENFIT acceptance was
# built with: matEffects OFF + SetBackExtrapToAxis + CATIMA vertex-gap eloss (rho 3.308e-5,
# matA = 1 for H2). The surviving fits in ~/a1954_C14_fit_300torr_matfx are matEffects=kTRUE
# with the RT geometry, a different configuration -- correcting those with this acceptance
# would mismatch data and simulation.
# The corrected energy lands in GetKinematicsXtr(), so the Ex cache is built with useXtr=kTRUE.
set -eo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); UKF=$(cd "$HERE/.." && pwd)
REPO=$(cd "$UKF/../../../.." && pwd)
IN=/home/yassid/a1954_C14_fit_300torr/in; OUT=/mnt/f/a1954_C14_gf_xtr; mkdir -p "$OUT"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
R=$1
[ -f "$OUT/$R.marker" ] && { echo "$R already done"; exit 0; }
ln -sf "$IN/${R}_reco.root" "$OUT/${R}_reco.root"
root -b -q -l "$UKF/pipeline/fitGenfit_C14.C(\"$R\",-1,\"$OUT/\",\"\",\"$OUT/\",-2.85,2,5,\"\",4.0,10.0,170.0,kFALSE,kFALSE,\"proton\",\"ATTPC_H300torr\",kTRUE,kTRUE,3.308e-5,1)" > "$OUT/${R}_fit.log" 2>&1
[ -s "$OUT/${R}_genfit.root" ] || { echo "$R FIT_FAILED"; exit 1; }
grep -q "CATIMA" "$OUT/${R}_fit.log" || { echo "$R CATIMA_NOT_ENABLED"; exit 1; }
echo COMPLETED > "$OUT/$R.marker"; echo "$R done"
