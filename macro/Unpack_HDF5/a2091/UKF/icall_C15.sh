#!/bin/bash
# Extract the ION CHAMBER for every a2091 15C run, both HDF5 formats.
#
# Per run: unpack the FRIB group to a TEMP <run>_FRIB.root, reduce it to the tiny
# <run>_ic.root summary (icmax + npulse per event), then DELETE the temp. Keeping the
# full FRIB output would be ~54 GB over 39 runs (8 channels x 2048 samples/event) while
# the analysis uses only two numbers per event -- so peak disk here is a few GB and the
# permanent product is ~25 MB.
#
# Format is auto-detected by unpackFRIB_C15.C:
#   legacy remerged (28 runs) : /frib/evt/evt<N>_1903          -> AtFRIBHDFUnpacker
#   raw merger      (11 runs) : /events/event_<N>/frib_physics/1903 -> AtMergerFRIBHDFUnpacker
#
#   ./icall_C15.sh              # all runs in runs_pp.txt, 4 parallel
#   ./icall_C15.sh "run_0177" 1
set -u
REPO=/home/yassid/fair_install/ATTPCROOTv2
HERE=$REPO/macro/Unpack_HDF5/a2091/UKF
IN="/media/yassid/Seagate Hub/ATTPC/Data/a2091/"
OUT=/home/yassid/a2091_C15_ic
TMP=$OUT/tmp
LOG=$OUT/logs
NPAR="${2:-4}"
mkdir -p "$TMP" "$LOG"

trap 'pkill -P $$ 2>/dev/null' EXIT INT TERM   # kill the xargs pool with the driver

set +u
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
set -u

# run_0140 is a corrupted HDF5 (bad symbol-table node -> core dump); skip it explicitly.
RUNS="${1:-$(tr ' ' '\n' < "$HERE/runs_pp.txt" | grep -v '^$' | grep -v run_0140 | tr '\n' ' ')}"
echo "[$(date +%H:%M:%S)] IC extraction over $(echo $RUNS | wc -w) runs, NPAR=$NPAR"

one(){ local r="$1"
  [ -f "$OUT/${r}_ic.root" ] && { echo "skip $r (done)"; return; }
  [ -f "$IN/${r}.h5" ]      || { echo "MISSING input $r"; return; }
  root -b -q -l "$HERE/pipeline/unpackFRIB_C15.C(\"$r\",\"$TMP/\")" > "$LOG/${r}_frib.log" 2>&1
  if [ ! -s "$TMP/${r}_FRIB.root" ]; then
    echo "[$(date +%H:%M:%S)] $r FRIB FAILED (see $LOG/${r}_frib.log)"; rm -f "$TMP/${r}_FRIB.root"; return
  fi
  root -b -q -l "$HERE/pipeline/icsum_C15.C(\"$r\",\"$TMP/\",\"$OUT/\")" > "$LOG/${r}_ic.log" 2>&1
  local line; line=$(grep -E "^${r} " "$LOG/${r}_ic.log" | head -1)
  rm -f "$TMP/${r}_FRIB.root"          # the 36 kB/event traces are not needed again
  if [ -f "$OUT/${r}_ic.root" ]; then echo "[$(date +%H:%M:%S)] ${line:-$r ok}"
  else echo "[$(date +%H:%M:%S)] $r SUMMARY FAILED"; fi
}
export -f one; export HERE IN OUT TMP LOG
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}

echo "[$(date +%H:%M:%S)] IC DONE: $(ls "$OUT"/*_ic.root 2>/dev/null | wc -l) summaries, $(du -sh "$OUT" 2>/dev/null | cut -f1) total"
rmdir "$TMP" 2>/dev/null
