#!/bin/bash
# GENFIT pass over the ALREADY-GATED (p,d) inputs written by fitpipe_pd_Be12.sh, so the
# browser explorer can switch between UKF and GENFIT on identical events/tracks.
# thetaMin is lowered to 5 deg (default 10) because the (p,d) deuterons are FORWARD.
# matEffects=kTRUE since the AtGenfitter dE/dx fix (5f87e625); geometry defaults to
# ATTPC_H300torr_RT (3.308e-5 g/cm3) as of 2026-08-25. NOTE: every (p,d) genfit set produced
# before that date ran matEffects=kTRUE against ATTPC_H600torr, i.e. 2x the real material,
# and unlike the (p,p') matFX-off sets that error was LIVE. Those results must be redone.
# Output dir bumped so the pre-fix results stay available.
#   ./genfit_pd_Be12.sh "run_0143 run_0147" 4
RUNS="${1:-run_0143}"; NPAR="${2:-4}"
REPO="/home/yassid/fair_install/ATTPCROOTv2-OpenKF"
HERE="$REPO/macro/Unpack_HDF5/a1954_Be12/UKF"
FITDIR="/home/yassid/a1954_Be12_fit_pd/"; IN="${FITDIR}in/"; LOG="${FITDIR}logs"
OUTDIR="${OUTDIR:-/home/yassid/a1954_Be12_genfix_pd/}"; mkdir -p "$OUTDIR"

source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"

one(){ local r="$1"; local L="$LOG/${r}_genfit.log"
  [ -f "$IN/${r}_reco.root" ] || { echo "no gated input $r"; return; }
  root -b -q -l "$HERE/pipeline/fitGenfit_Be12.C(\"$r\",-1,\"$IN\",\"\",\"$OUTDIR\",-2.85,2,5,\"\",4.0,5.0,170.0,kTRUE,kFALSE,\"deuteron\")" > "$L" 2>&1
  echo "[$(date +%H:%M:%S)] $r genfit=$([ -f $OUTDIR${r}_genfit.root ]&&echo ok||echo FAIL)"
}
export -f one; export HERE FITDIR IN LOG OUTDIR
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "[$(date +%H:%M:%S)] (p,d) GENFIT PASS DONE -> $FITDIR"
