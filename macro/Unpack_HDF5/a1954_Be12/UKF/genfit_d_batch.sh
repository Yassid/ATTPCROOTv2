#!/bin/bash
RUNS="$1"; NPAR="${2:-4}"; DIR="/home/yassid/a1954_Be12_reco_hdb"
HERE="$(cd "$(dirname "$0")" && pwd)"; LOG="$DIR/logs"; mkdir -p "$LOG"
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh
source /home/yassid/fair_install/ATTPCROOTv2/build_genfit/config.sh >/dev/null 2>&1
one(){ local r="$1"
  [ -f "$DIR/${r}_reco.root" ] || { echo "skip $r"; return; }
  [ -f "$DIR/${r}_genfit_d.root" ] && { echo "skip $r (exists)"; return; }
  # fitGenfit_Be12(file,nEv,ioDir,outSuffix,outDir,bField,minIter,maxIter,pidGate,measSigma,thMin,thMax,matEff,backSeed,particle)
  root -b -q -l "$HERE/pipeline/fitGenfit_Be12.C(\"$r\",-1,\"$DIR/\",\"_d\",\"$DIR/\",-2.85,2,5,\"\",4.0,10,170,kTRUE,kFALSE,\"deuteron\")" > "$LOG/${r}_genfitd.log" 2>&1
  echo "[$(date +%H:%M:%S)] genfit_d $r exit $?"; }
export -f one; export DIR HERE LOG
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "GENFIT deuteron batch done."
