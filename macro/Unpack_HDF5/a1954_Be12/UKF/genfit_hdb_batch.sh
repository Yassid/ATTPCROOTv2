#!/bin/bash
# GENFIT fit over HDBSCAN reco -> <run>_genfit.root (proton). Uses build_genfit/.
RUNS="$1"; NPAR="${2:-4}"; DIR="/home/yassid/a1954_Be12_reco_hdb"
HERE="$(cd "$(dirname "$0")" && pwd)"; LOG="$DIR/logs"; mkdir -p "$LOG"
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh
source /home/yassid/fair_install/ATTPCROOTv2/build_genfit/config.sh >/dev/null 2>&1
one(){ local r="$1"
  [ -f "$DIR/${r}_reco.root" ] || { echo "skip $r (no reco)"; return; }
  [ -f "$DIR/${r}_genfit.root" ] && { echo "skip $r (genfit exists)"; return; }
  root -b -q -l "$HERE/pipeline/fitGenfit_Be12.C(\"$r\",-1,\"$DIR/\",\"\",\"$DIR/\")" > "$LOG/${r}_genfit.log" 2>&1
  echo "[$(date +%H:%M:%S)] genfit $r exit $?"; }
export -f one; export DIR HERE LOG
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "GENFIT hdb batch done."
