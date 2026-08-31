#!/bin/bash
# Bulk-unpack the Dec2014 alpha runs with the FIXED fr19port chain (AsAd merge corrected,
# validated digit-for-digit against the 2014 reference analysis). Hits only: the decoder
# task does not persist AtRawEvent (~30 KB/event vs ~1.3 MB).
#
# The macro needs an exact event count -- asking for more events than the run holds CRASHES
# at end of data (rc=129) and truncates the output, so count first with count_events.C.
#   usage: unpack_all_fr19.sh <run> [run ...]
source /home/yassid/fair_install/ATTPCROOTv2_fr19port/setup_fr19port.sh > /dev/null 2>&1
cd /home/yassid/fair_install/ATTPCROOTv2_fr19port/macro/Unpack_GETDecoder2 || exit 1
OUT=/home/yassid/dec2014_alphas_reco
LOG=$OUT/logs
mkdir -p "$OUT" "$LOG"

for run in "$@"; do
  o="$OUT/alpha_${run}_hits.root"
  [ -s "$o" ] && { echo "[skip] $run already done"; continue; }
  list="runfiles/NSCL/Dec2014_alphas/alpha_${run}.txt"
  [ -s "$list" ] || { echo "[skip] $run : no runfile"; continue; }

  n=$(timeout 1800 root -l -b -q -e 'gSystem->Load("libAtReconstruction");' \
        "count_events.C(\"$list\")" 2>/dev/null | grep -oE "NEVENTS [0-9]+" | grep -oE "[0-9]+")
  if [ -z "$n" ] || [ "$n" -le 0 ]; then echo "[skip] $run : could not count events"; continue; fi

  echo "[run ] $run  events=$n  $(date +%H:%M:%S)"
  timeout 7200 root -l -b -q \
    "run_unpack_Dec2014_alphas.C(\"$list\",\"$o\",$n,\"Lookup20141208.xml\")" \
    > "$LOG/$run.log" 2>&1
  rc=$?
  sz=$(du -h "$o" 2>/dev/null | cut -f1)
  pads=$(grep -oE "Number of Pads : [0-9]+" "$LOG/$run.log" | tail -1 | grep -oE "[0-9]+$")
  echo "[done] $run rc=$rc events=$n pads/event=${pads:-?} size=${sz:-none}  $(date +%H:%M:%S)"
done
echo "ALL_DONE $(date +%H:%M:%S)"
