#!/bin/bash
# Build the slim LOCAL AtPatternEvent cache from the big reco.root files on F.
# One drvfs pass; afterwards all PID/fit passes read locally (~250x faster).
#   ./slim_cache_batch.sh "run_0138 run_0056 ..." 4
RUNS="${1:-run_0138}"; NPAR="${2:-4}"
REPO="/home/yassid/fair_install/ATTPCROOTv2"
HERE="$REPO/macro/Unpack_HDF5/a2091/UKF"
IN="/home/yassid/a2091_C15_reco/"
OUT="/home/yassid/a2091_C15_reco_slim/"; LOG="${OUT}logs"; mkdir -p "$LOG"

source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"

slim_one(){ local r="$1"
  [ -f "$OUT/${r}_slim.root" ] && { echo "skip $r"; return; }
  [ -f "$IN/${r}_reco.root" ] || { echo "no input $r"; return; }
  root -b -q -l "$HERE/pipeline/slim_cache_C15.C(\"$r\",\"$IN\",\"$OUT\")" > "$LOG/${r}_slim.log" 2>&1
  echo "[$(date +%H:%M:%S)] $r exit $? ($(ls -la $OUT/${r}_slim.root 2>/dev/null | awk '{printf "%d MB",$5/1e6}'))"
}
export -f slim_one; export OUT IN HERE LOG
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'slim_one "$@"' _ {}
echo "[$(date +%H:%M:%S)] slim cache DONE  ($(du -sh $OUT 2>/dev/null|cut -f1) total)"
