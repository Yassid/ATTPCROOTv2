#!/bin/bash
RUNS="${1}"; NPAR="${2:-1}"
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="/home/yassid/a1954_Be12_reco"; LOG="$OUT/logs"; mkdir -p "$LOG"
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh
source /home/yassid/fair_install/ATTPCROOTv2/build/config.sh >/dev/null 2>&1
fr(){ local r="$1"
  [ -f "$OUT/${r}_FRIB.root" ] && { echo "skip $r"; return; }
  root -b -q -l "$HERE/pipeline/unpackFRIB_Be12.C(\"$r\")" > "$LOG/${r}_frib.log" 2>&1
  echo "[$(date +%H:%M:%S)] FRIB $r exit $? ($(ls -la $OUT/${r}_FRIB.root 2>/dev/null|awk '{printf "%d MB",$5/1e6}'))"
}
export -f fr; export OUT HERE LOG
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'fr "$@"' _ {}
echo "FRIB batch done."
