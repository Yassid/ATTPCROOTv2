#!/bin/bash
# Bulk-unpack the curated Dec2014 alpha runs, hits only (no AtRawEvent -> ~32 MB/run).
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh
source /home/yassid/fair_install/ATTPCROOTv2/build/config.sh > /dev/null 2>&1
cd /home/yassid/fair_install/ATTPCROOTv2/macro/Unpack_GRAW_Dec2014 || exit 1
OUT=/home/yassid/dec2014_alphas_reco
LOG=$OUT/logs
mkdir -p "$OUT" "$LOG"

for run in "$@"; do
  o="$OUT/alpha_${run}_hits.root"
  if [ -s "$o" ]; then echo "[skip] $run already done"; continue; fi
  read -r listfile ncobo < <(/home/yassid/dec2014_links/make_links.sh "$run")
  if [ "$listfile" = "NO_DIR" ] || [ "$listfile" = "NO_DATA" ]; then
    echo "[skip] $run : $listfile"; continue
  fi
  echo "[run ] $run  cobos=$ncobo  $(date +%H:%M:%S)"
  timeout 3600 root -l -b -q \
    "unpack_Dec2014_alphas.C(\"$listfile\",\"$o\",0,$ncobo,\"Lookup20141208.xml\",\"ATTPC.alpha.par\",kFALSE)" \
    > "$LOG/$run.log" 2>&1
  rc=$?
  n=$(grep -oE "Unpacking [0-9]+ events" "$LOG/$run.log" | grep -oE "[0-9]+" | head -1)
  sz=$(du -h "$o" 2>/dev/null | cut -f1)
  echo "[done] $run rc=$rc events=${n:-?} size=${sz:-none}"
done
echo "ALL_DONE $(date +%H:%M:%S)"
