#!/usr/bin/env bash
# Drive the PID accumulation over several seeds, N at a time, under a memory tracer with a guard.
#
#   ./run_pid_accumulation.sh [-j N] [seed ...]        default: -j 4, seeds 2001..2006
#
# WHY THE GUARD EXISTS. Two earlier attempts at this job took the whole WSL VM down, both during
# reco with three seeds running concurrently, and neither left any evidence: when the VM dies the
# dmesg goes with it and all that remains is a truncated ROOT file. Nothing was watching memory.
# So the tracer samples to a file on F: (which survives the VM) and kills the largest ROOT child
# at 85 % of MemTotal rather than letting the host get there.
#
# WHY 4 IS A REASONABLE DEFAULT, AND WHAT WOULD CHANGE IT. A single reco measured ~1.0 GB at
# event 150 rising to ~1.15 GB at event 2600, i.e. ~0.07 MB/event, which projects to under 2 GB
# over the full 12000. Four of those is ~8 GB against a 23.5 GB cap on 4 of 8 processors. NOTE
# that this also means three concurrent recos come to ~6 GB and do NOT by themselves explain the
# crashes -- so the cause is not established, and the honest reading is that either a later spike
# exists that the short profile never reached, or something outside these jobs was involved.
# Raise -j only after a full-length profile shows no spike; memtrace.csv now records the SUM over
# jobs precisely so that question can be answered from the run itself.
#
# A failed seed does not stop the others: accumulate_pid.sh is resumable per stage on evidence
# that a stage finished, so a seed picks up where it stopped on the next invocation.
set -uo pipefail
JOBS=4
if [ "${1:-}" = "-j" ]; then JOBS=${2:?-j needs a number}; shift 2; fi
SEEDS=("$@"); [ ${#SEEDS[@]} -eq 0 ] && SEEDS=(2001 2002 2003 2004 2005 2006)
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUT=/mnt/f/a1975_C16_pp_pid
NEV=12000
TRACE="$OUT/memtrace.csv"
STAGEFILE="$OUT/.stage"
LOG="$OUT/accum_run.log"

mkdir -p "$OUT"
echo "j$JOBS" > "$STAGEFILE"
# 2 s rather than 5: with several jobs allocating at once the guard has less slack, and the
# sampling cost is nil next to the reco.
"$DIR/memtrace.sh" "$TRACE" "$STAGEFILE" 85 2 &
TRACER=$!

# A bare `trap ... TERM` that only cleans up does NOT stop the run: bash executes the handler and
# then carries straight on with the next seed, which is exactly what happened the first time this
# was stopped. The signal handler must therefore tear down the children and exit. Killing by
# pattern is safe here because neither pattern matches this script's own command line.
cleanup() { kill "$TRACER" 2>/dev/null; echo idle > "$STAGEFILE"; }
stop() {
   cleanup
   # Kill the whole job tree, not just xargs: xargs dying does not touch the ROOT process its
   # child already started, and that ROOT process is the one holding the memory.
   [ -n "${XARGS:-}" ] && kill -9 "$XARGS" 2>/dev/null
   pkill -9 -f 'accumulate_pid\.sh' 2>/dev/null
   pkill -9 -f 'root\.exe.*(a1975_C16_pp_pid|C16_pp_a1975_sim)' 2>/dev/null
   echo "[$(date +%F' '%H:%M:%S)] === stopped by signal ===" | tee -a "$LOG"
   exit 130
}
trap cleanup EXIT
trap stop INT TERM

say() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$LOG"; }
say "=== accumulation start: -j $JOBS, seeds ${SEEDS[*]} (tracer $TRACER, guard 85%) ==="

# xargs runs in the BACKGROUND and is waited on, rather than in the foreground. bash defers a
# trap until the current foreground child returns, so with xargs in front a SIGTERM to this
# script did nothing until the whole run was over -- which is exactly what happened when the run
# was first stopped: the handler never ran, and the driver calmly started the next seed. With
# wait, the signal is handled at once.
printf '%s\n' "${SEEDS[@]}" | xargs -P "$JOBS" -I{} \
   bash -c 'start=$SECONDS
      if "$0" "$1" "$2" "$3" >> "$4" 2>&1; then
         echo "[$(date +%F" "%H:%M:%S)] seed $1 COMPLETED in $(( (SECONDS-start)/60 )) min"
      else
         echo "[$(date +%F" "%H:%M:%S)] seed $1 FAILED after $(( (SECONDS-start)/60 )) min (resumable)"
      fi | tee -a "$4"' \
   "$DIR/accumulate_pid.sh" {} "$NEV" "$OUT" "$LOG" &
XARGS=$!
wait "$XARGS"

done_n=$(ls "$OUT"/*.marker 2>/dev/null | wc -l)
say "=== finished: $done_n of ${#SEEDS[@]} seeds carry a COMPLETED marker ==="
say "=== peak: $(awk -F, 'NR>1 && $6+0>m {m=$6; l=$0} END{print m"% used -- "l}' "$TRACE") ==="
