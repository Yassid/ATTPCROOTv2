#!/bin/bash
# Unpack + PSA for one Dec2014 alpha run. Hits only (pulses not persisted).
#   usage: unpack_one.sh <parFile> <run>
PAR=$1; run=$2
source /home/yassid/fair_install/ATTPCROOTv2_fr19port/setup_fr19port.sh > /dev/null 2>&1
cd /home/yassid/fair_install/ATTPCROOTv2_fr19port/macro/Unpack_GETDecoder2 || exit 1
OUT=/home/yassid/dec2014_alphas_reco/lowP; LOG=$OUT/logs
mkdir -p "$OUT" "$LOG"
rr=$(printf "run_%04d" "$run")
o="$OUT/alpha_${rr}_hits.root"
[ -s "$o" ] && { echo "[skip] $rr already done"; exit 0; }
list="runfiles/NSCL/Dec2014_alphas/alpha_${rr}.txt"
[ -s "$list" ] || { echo "[skip] $rr : no runfile"; exit 0; }

# The macros print per-frame chatter for every event; unfiltered that reached 39 GB for a
# single run. Keep only lines that matter.
FILTER='fCoboFrameInfo|fFrameInfoIdx|Returned event ID|Not full in|fTargetFrameInfoIdx|Event Number|Valid pads|MC Simulated'

n=$(timeout 3600 root -l -b -q -e 'gSystem->Load("libAtReconstruction");' \
      "count_events.C(\"$list\")" 2>/dev/null | grep -oE "NEVENTS [0-9]+" | grep -oE "[0-9]+")
if [ -z "$n" ] || [ "$n" -le 0 ]; then echo "[skip] $rr : could not count events"; exit 1; fi

echo "[run ] $rr events=$n $(date +%H:%M:%S)"
timeout 14400 root -l -b -q \
  "run_unpack_Dec2014_alphas.C(\"$list\",\"$o\",$n,\"Lookup20141208.xml\",\"$PAR\")" 2>&1 \
  | grep -avE "$FILTER" | tail -400 > "$LOG/$rr.log"
rc=${PIPESTATUS[0]}
sz=$(du -h "$o" 2>/dev/null | cut -f1)
echo "[done] $rr rc=$rc events=$n size=${sz:-none} $(date +%H:%M:%S)"
