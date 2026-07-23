#!/bin/bash
# UKF fit over HDBSCAN reco -> <run>_ukf.root (proton hypothesis). Uses default build/.
RUNS="$1"; NPAR="${2:-6}"; PARTS="${3:-proton}"; DIR="/home/yassid/a1954_Be12_reco_hdb"
HERE="$(cd "$(dirname "$0")" && pwd)"; LOG="$DIR/logs"; mkdir -p "$LOG"
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh
source /home/yassid/fair_install/ATTPCROOTv2/build/config.sh >/dev/null 2>&1
one(){ local r="$1"
  [ -f "$DIR/${r}_reco.root" ] || { echo "skip $r (no reco)"; return; }
  [ -f "$DIR/${r}_ukf.root" ] && { echo "skip $r (ukf exists)"; return; }
  root -b -q -l "$HERE/pipeline/fitUKF_Be12.C(\"$r\",-1,\"$PARTS\",-1,2.85,6.5e-5,\"\",\"$DIR/\")" > "$LOG/${r}_ukf.log" 2>&1
  echo "[$(date +%H:%M:%S)] ukf $r exit $?"; }
export -f one; export DIR HERE PARTS LOG
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "UKF hdb batch done."
