#!/bin/bash
RUNS="$1"; NPAR="${2:-6}"; DIR="/home/yassid/a1954_Be12_reco_hdb"
HERE="$(cd "$(dirname "$0")" && pwd)"; LOG="$DIR/logs"; mkdir -p "$LOG"
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh
source /home/yassid/fair_install/ATTPCROOTv2/build/config.sh >/dev/null 2>&1
one(){ local r="$1"
  [ -f "$DIR/${r}_reco.root" ] || { echo "skip $r"; return; }
  [ -f "$DIR/${r}_ukf_d.root" ] && { echo "skip $r (exists)"; return; }
  # fitUKF_Be12(file,nEv,particles,bSign,bMag,gasDens,outSuffix,ioDir)
  root -b -q -l "$HERE/pipeline/fitUKF_Be12.C(\"$r\",-1,\"deuteron\",-1,2.85,6.5e-5,\"_d\",\"$DIR/\")" > "$LOG/${r}_ukfd.log" 2>&1
  echo "[$(date +%H:%M:%S)] ukf_d $r exit $?"; }
export -f one; export DIR HERE LOG
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "UKF deuteron batch done."
