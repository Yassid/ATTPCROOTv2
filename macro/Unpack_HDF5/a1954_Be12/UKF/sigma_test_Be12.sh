#!/bin/bash
# Covariance test: refit the SAME gated tracks with the sigmas swapped (UKF 0.5->4.0 mm,
# GENFIT 4.0->0.5 mm) to separate chi2-scale effects from real fitter differences.
# See pp/FITTER_COVARIANCE_TEST.md for the results.
# matched-sigma test: UKF re-run at 4.0 mm, GENFIT re-run at 0.5 mm, same gated tracks
REPO="/home/yassid/fair_install/ATTPCROOTv2-OpenKF"
HERE="$REPO/macro/Unpack_HDF5/a1954_Be12/UKF"
IN="/home/yassid/a1954_Be12_fit/in/"; OUT="/home/yassid/a1954_Be12_sigtest/"; mkdir -p "$OUT/logs"
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
one(){ r="$1"
  root -b -q -l "$HERE/pipeline/fitUKF_Be12.C(\"$r\",-1,\"proton\",-1,2.85,6.5e-5,\"_s4\",\"$IN\",4.0,0.1,1,10,\"$OUT\")" > "$OUT/logs/${r}_u4.log" 2>&1
  root -b -q -l "$HERE/pipeline/fitGenfit_Be12.C(\"$r\",-1,\"$IN\",\"_s05\",\"$OUT\",-2.85,2,5,\"\",0.5,10.0,170.0,kFALSE,kFALSE,\"proton\")" > "$OUT/logs/${r}_g05.log" 2>&1
  echo "$r  ukf_s4=$([ -f $OUT${r}_ukf_s4.root ]&&echo ok||echo FAIL)  genfit_s05=$([ -f $OUT${r}_genfit_s05.root ]&&echo ok||echo FAIL)"
}
export -f one; export HERE IN OUT
printf "%s\n" $1 | xargs -P 4 -I{} bash -c 'one "$@"' _ {}
