#!/usr/bin/env bash
# Drive the 46Ar(3He,d)47K samples, N at a time, under a memory tracer with a guard.
#
#   ./run_3Hed_accumulation.sh [-j N] [state:seed ...]
#
# Default job list is the three 47K levels with two seeds each:
#     gs:3001 gs:3002   360:3011 360:3012   2020:3021 2020:3022
# Seeds are distinct ACROSS states as well as within them, so no two samples anywhere in the set
# share a random sequence.
#
# WHY THE GUARD EXISTS. Two a1975 runs took the whole WSL VM down during reco, and when the VM
# dies the dmesg goes with it -- all that survives is a truncated ROOT file. So the tracer samples
# to a file on F: (which survives the VM) and kills the largest ROOT child at 85 % of MemTotal
# rather than letting the host get there. The tracer itself is a1975's, reused rather than copied
# so there is one implementation of the guard to trust.
#
# WHY -j 3. Measured on this chain, not assumed: a single reco sat at 1.14 GB after 140 events of
# a 400-event sample. Three of those is ~3.5 GB against a 23.5 GB cap. That is deliberately more
# conservative than a1975's -j 4, for a reason specific to this simulation: the 46Ar beam is
# Z = 18 and deposits about nine times the charge per unit length that a1975's 16C did, so the
# events are heavier -- ~1e6 points per event here -- and the memory profile over a full 12000
# entry sample has not been measured. Raise it only after watching memtrace.csv through a whole
# sample.
#
# COST, measured: reco runs ~0.46 s/entry, so a 12000-entry sample (6000 reactions) is about
# 90 minutes; generation is ~3 minutes and the PID pass about a minute. Six samples at -j 3 is
# therefore roughly 3 hours.
#
# A failed sample does not stop the others: accumulate_3Hed.sh is resumable per stage, so it picks
# up where it stopped on the next invocation.
set -uo pipefail
JOBS=3
if [ "${1:-}" = "-j" ]; then JOBS=${2:?-j needs a number}; shift 2; fi
JOBLIST=("$@")
[ ${#JOBLIST[@]} -eq 0 ] && JOBLIST=(gs:3001 gs:3002 360:3011 360:3012 2020:3021 2020:3022)
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
# Configuration comes from the environment so one driver serves the whole field/pad matrix:
#   OUT     where reco+pid land          BT   field in tesla
#   SIMDIR  where the sims live/go       PAD  pad pitch in mm, <=0 = real AT-TPC plane
# Generation depends on the FIELD ONLY, so the 2 mm configurations point SIMDIR at the sims the
# corresponding field already produced instead of regenerating them.
OUT=${OUT:-/mnt/f/ar46_3hed}
BT=${BT:-2.85}
PAD=${PAD:--1}
SIMDIR=${SIMDIR:-$OUT}
NEV=${NEV:-12000}
TRACE="$OUT/memtrace.csv"
STAGEFILE="$OUT/.stage"
LOG="$OUT/accum_run.log"

mkdir -p "$OUT"
echo "j$JOBS" > "$STAGEFILE"
"$DIR/../16C_pp_a1975/memtrace.sh" "$TRACE" "$STAGEFILE" 85 2 &
TRACER=$!

# A bare `trap ... TERM` that only cleans up does NOT stop the run: bash runs the handler and then
# carries straight on with the next job. The handler must tear down the children and exit.
cleanup() { kill "$TRACER" 2>/dev/null; echo idle > "$STAGEFILE"; }
stop() {
   cleanup
   # Kill the whole job tree, not just xargs: xargs dying does not touch the ROOT process its
   # child already started, and that is the process holding the memory.
   [ -n "${XARGS:-}" ] && kill -9 "$XARGS" 2>/dev/null
   pkill -9 -f 'accumulate_3Hed\.sh' 2>/dev/null
   pkill -9 -f 'root\.exe.*(ar46_3hed|Ar46_3Hed_sim|run_reco_Ar46)' 2>/dev/null
   echo "[$(date +%F' '%H:%M:%S)] === stopped by signal ===" | tee -a "$LOG"
   exit 130
}
trap cleanup EXIT
trap stop INT TERM

say() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$LOG"; }
say "=== 3He,d accumulation start: -j $JOBS, B = $BT T, pads = $PAD mm, sims from $SIMDIR, out $OUT ==="

# xargs runs in the BACKGROUND and is waited on: bash defers a trap until the current foreground
# child returns, so with xargs in front a SIGTERM to this script does nothing until the whole run
# is over. With wait, the signal is handled at once.
printf '%s\n' "${JOBLIST[@]}" | xargs -P "$JOBS" -I{} \
   bash -c 'start=$SECONDS; state=${1%%:*}; seed=${1##*:}
      if "$0" "$state" "$seed" "$2" "$3" "$5" "$6" "$7" >> "$4" 2>&1; then
         echo "[$(date +%F" "%H:%M:%S)] $1 COMPLETED in $(( (SECONDS-start)/60 )) min"
      else
         echo "[$(date +%F" "%H:%M:%S)] $1 FAILED after $(( (SECONDS-start)/60 )) min (resumable)"
      fi | tee -a "$4"' \
   "$DIR/accumulate_3Hed.sh" {} "$NEV" "$OUT" "$LOG" "$BT" "$PAD" "$SIMDIR" &
XARGS=$!
wait "$XARGS"

done_n=$(ls "$OUT"/*.marker 2>/dev/null | wc -l)
say "=== finished: $done_n of ${#JOBLIST[@]} samples carry a COMPLETED marker ==="
say "=== peak: $(awk -F, 'NR>1 && $6+0>m {m=$6; l=$0} END{print m"% used -- "l}' "$TRACE") ==="
