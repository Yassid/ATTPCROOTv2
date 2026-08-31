#!/bin/bash
# Unpack + PSA only for the low-pressure (149-150 torr) Dec 2014 alpha runs.
# Hits only -- the decoder task does not persist AtRawEvent, so no pulses are saved.
#   usage: unpack_lowP.sh <parFile> <run> [run ...]
PAR=$1; shift
source /home/yassid/fair_install/ATTPCROOTv2_fr19port/setup_fr19port.sh > /dev/null 2>&1
cd /home/yassid/fair_install/ATTPCROOTv2_fr19port/macro/Unpack_GETDecoder2 || exit 1
OUT=/home/yassid/dec2014_alphas_reco/lowP
LOG=$OUT/logs
mkdir -p "$OUT" "$LOG"

for run in "$@"; do
  rr=$(printf "run_%04d" "$run")
  o="$OUT/alpha_${rr}_hits.root"
  [ -s "$o" ] && { echo "[skip] $rr already done"; continue; }
  list="runfiles/NSCL/Dec2014_alphas/alpha_${rr}.txt"
  [ -s "$list" ] || { echo "[skip] $rr : no runfile"; continue; }

  # Exact event count: overrunning end-of-data crashes and truncates the output.
  n=$(timeout 3600 root -l -b -q -e 'gSystem->Load("libAtReconstruction");' \
        "count_events.C(\"$list\")" 2>/dev/null | grep -oE "NEVENTS [0-9]+" | grep -oE "[0-9]+")
  if [ -z "$n" ] || [ "$n" -le 0 ]; then echo "[skip] $rr : could not count events"; continue; fi

  echo "[run ] $rr  events=$n  $(date +%H:%M:%S)"
  timeout 10800 root -l -b -q \
    "run_unpack_Dec2014_alphas.C(\"$list\",\"$o\",$n,\"Lookup20141208.xml\",\"$PAR\")" \
    > "$LOG/$rr.log" 2>&1
  rc=$?
  sz=$(du -h "$o" 2>/dev/null | cut -f1)
  pads=$(grep -oE "Number of Pads : [0-9]+" "$LOG/$rr.log" | tail -1 | grep -oE "[0-9]+$")
  echo "[done] $rr rc=$rc events=$n pads/event=${pads:-?} size=${sz:-none}  $(date +%H:%M:%S)"
done
echo "ALL_DONE $(date +%H:%M:%S)"
