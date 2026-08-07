#!/bin/bash
# Unpack the magnet-OFF calibration runs (96-113, 299.5 torr) -- hits only.
#   usage: unpack_noB.sh <run> [run ...]
PAR=ATTPC.alpha_300torr_noB.par
source /home/yassid/fair_install/ATTPCROOTv2_fr19port/setup_fr19port.sh > /dev/null 2>&1
cd /home/yassid/fair_install/ATTPCROOTv2_fr19port/macro/Unpack_GETDecoder2 || exit 1
OUT=/home/yassid/dec2014_alphas_reco/noB; LOG=$OUT/logs
mkdir -p "$OUT" "$LOG"
for run in "$@"; do
  rr=$(printf "run_%04d" "$run"); o="$OUT/alpha_${rr}_hits.root"
  [ -s "$o" ] && { echo "[skip] $rr already done"; continue; }
  list="runfiles/NSCL/Dec2014_alphas/alpha_${rr}.txt"
  [ -s "$list" ] || { echo "[skip] $rr : no runfile"; continue; }
  n=$(timeout 3600 root -l -b -q -e 'gSystem->Load("libAtReconstruction");' \
        "count_events.C(\"$list\")" 2>/dev/null | grep -oE "NEVENTS [0-9]+" | grep -oE "[0-9]+")
  [ -z "$n" ] || [ "$n" -le 0 ] && { echo "[skip] $rr : could not count events"; continue; }
  echo "[run ] $rr  events=$n  $(date +%H:%M:%S)"
  timeout 10800 root -l -b -q \
    "run_unpack_Dec2014_alphas.C(\"$list\",\"$o\",$n,\"Lookup20141208.xml\",\"$PAR\")" > "$LOG/$rr.log" 2>&1
  echo "[done] $rr rc=$? events=$n size=$(du -h "$o" 2>/dev/null|cut -f1)  $(date +%H:%M:%S)"
done
echo "ALL_DONE $(date +%H:%M:%S)"
