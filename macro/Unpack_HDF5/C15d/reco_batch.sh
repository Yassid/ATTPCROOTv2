#!/usr/bin/env bash
# C15d: reco + gain-matched PID cache, over a run list. Resumable, disk-guarded.
#
#   ./reco_batch.sh [nparallel] [runlist] [nEvents]
#
# Each run does ONE pass over the raw HDF5: unpack -> PSA -> SC -> PRA -> PID -> gain match,
# writing <run>_reco.root, then dumps the small <run>_pid.root plane from it. Both steps are
# skipped if their output already exists, so re-running picks up where it stopped.
#
# NPARALLEL: 8 on this machine (32 cores). MEASURED IN EVENTS PER SECOND:
#
#     NPAR=3    49 events/s aggregate
#     NPAR=8   162 events/s aggregate
#
# ★ MEASURE THIS IN EVENTS/S, NOT MB/S. Reading MB/s out of /proc/diskstats is worse than
# useless here: with ~50 GB of page cache holding recently-read parts of the raw files, sda can
# show 1 MB/s while every job is happily processing thousands of events. Judged on that metric
# this stage looks I/O bound and looks like it degrades with parallelism -- both conclusions are
# artifacts of the cache. Judged on events/s it scales nearly linearly from 3 to 8.
#
# Some jobs still sit in D state and a slow one can be starved for minutes while others stream;
# that is normal here and it does not mean the pool is too wide -- the starved run catches up
# when a neighbour finishes. Total for the 105-run set at 162 ev/s is roughly 5 h.
#
# Traps this script exists to avoid, all of which have bitten this workspace before:
#   - sourcing config.sh under `set -u` kills the shell silently -> `set +u` around it
#   - `ps -eo args | grep <pat>` and `pkill -f <pat>` MATCH THEMSELVES -> pgrep -f / bracket
#   - killing the driver does NOT kill an xargs pool -> trap pkill -P $$
#   - a half-written reco is happily "recovered" by ROOT and reports plausible numbers ->
#     write to <run>_reco.root.part and mv only on success, so a partial file is never
#     mistaken for a finished one

set -eo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"

NPAR="${1:-16}"
# runs_d2.txt, NOT runs_d.txt. runs_d.txt is the original list and still carries the 28 runs
# >=106, which were later established to be a HYDROGEN target (dE/dx ratio 0.52 against 0.500
# predicted for H2/D2, and arclength jumping 162-166 -> 182-185 mm). runs_d2.txt is the curated
# D2-only set written once that was known. Leaving the default on the stale list reconstructed
# all 28 hydrogen runs, and because they are the LARGEST files a downstream `ls -S | head` then
# fitted 7 of them as deuterium -- which put the (d,d) beam energy at 54 MeV instead of 190.
RUNLIST="${2:-$HERE/runs_d2.txt}"
NEVENTS="${3:--1}"

RECO_DIR="${C15D_RECO:-/home/yassid/C15d_reco}"
LOG_DIR="${C15D_LOGS:-/home/yassid/C15d_logs}"
MIN_FREE_GB="${MIN_FREE_GB:-40}"

set +u
# shellcheck disable=SC1091
source "$REPO/build/config.sh" >/dev/null 2>&1
set -u

if [[ -z "${VMCWORKDIR:-}" ]]; then
   echo "ERROR: config.sh did not set VMCWORKDIR ($REPO/build/config.sh)" >&2
   exit 1
fi
if [[ ! -f "$RUNLIST" ]]; then
   echo "ERROR: run list not found: $RUNLIST" >&2
   exit 1
fi

mkdir -p "$RECO_DIR" "$LOG_DIR"

mapfile -t RUNS < <(grep -vE '^\s*(#|$)' "$RUNLIST")
if [[ ${#RUNS[@]} -eq 0 ]]; then
   echo "ERROR: run list $RUNLIST is empty. Refusing to run -- an empty list here has" >&2
   echo "       previously turned the whole stage into a silent no-op." >&2
   exit 1
fi

echo "=== C15d reco_batch ==="
echo "  runs      : ${#RUNS[@]} from $RUNLIST"
echo "  parallel  : $NPAR"
echo "  nEvents   : $NEVENTS  (-1 = all)"
echo "  reco dir  : $RECO_DIR"
echo "  free      : $(df -BG --output=avail "$RECO_DIR" | tail -1 | tr -d ' G') GB (guard ${MIN_FREE_GB} GB)"
echo

# Kill the whole process group on exit: xargs keeps spawning children after the parent dies.
trap 'pkill -P $$ 2>/dev/null || true' EXIT INT TERM

do_run() {
   local run="$1"
   local reco="$RECO_DIR/${run}_reco.root"
   local pid="$RECO_DIR/${run}_pid.root"
   local log="$LOG_DIR/${run}.log"

   local free_gb
   free_gb=$(df -BG --output=avail "$RECO_DIR" | tail -1 | tr -d ' G')
   if (( free_gb < MIN_FREE_GB )); then
      echo "[$run] SKIP: only ${free_gb} GB free (< ${MIN_FREE_GB})" | tee -a "$log"
      return 0
   fi

   if [[ -s "$reco" ]]; then
      echo "[$run] reco exists, skipping"
   else
      echo "[$run] reco ..."
      # Write into a staging dir and mv into place ONLY on success. A file size test is not
      # enough on its own: ROOT happily "recovers" a truncated reco and reports plausible
      # numbers, so a job killed mid-write would leave a non-empty file that the skip test
      # above accepts as finished, and the run silently enters the analysis short. mv within
      # the same filesystem is atomic, so a reco under its final name is always complete.
      local part_dir="$RECO_DIR/.part/$run"
      rm -rf "$part_dir"
      mkdir -p "$part_dir"
      if root -b -q "$HERE/unpackReco_C15d.C(\"$run\", $NEVENTS, false, \"$part_dir/\")" \
            >"$log" 2>&1 && [[ -s "$part_dir/${run}_reco.root" ]]; then
         mv -f "$part_dir/${run}_reco.root" "$reco"
         rm -rf "$part_dir"
         echo "[$run] reco OK  ($(du -h "$reco" | cut -f1))"
      else
         echo "[$run] RECO FAILED -- see $log"
         rm -rf "$part_dir"
         return 0
      fi
   fi

   if [[ -s "$pid" ]]; then
      echo "[$run] pid exists, skipping"
   else
      local pid_part="$RECO_DIR/.part/${run}_pidcache"
      rm -rf "$pid_part"
      mkdir -p "$pid_part"
      if root -b -q "$HERE/pidntuple_C15d.C(\"$run\", \"$RECO_DIR/\", \"$pid_part/\")" >>"$log" 2>&1 \
            && [[ -s "$pid_part/${run}_pid.root" ]]; then
         mv -f "$pid_part/${run}_pid.root" "$pid"
         rm -rf "$pid_part"
         echo "[$run] pid OK   ($(du -h "$pid" | cut -f1))"
      else
         echo "[$run] PID CACHE FAILED -- see $log"
         rm -rf "$pid_part"
      fi
   fi
}
export -f do_run
export HERE RECO_DIR LOG_DIR MIN_FREE_GB NEVENTS

printf '%s\n' "${RUNS[@]}" | xargs -P "$NPAR" -I{} bash -c 'do_run "$@"' _ {}

echo
echo "=== done: $(ls -1 "$RECO_DIR"/*_reco.root 2>/dev/null | wc -l) recos, "\
"$(ls -1 "$RECO_DIR"/*_pid.root 2>/dev/null | wc -l) PID caches in $RECO_DIR ==="
