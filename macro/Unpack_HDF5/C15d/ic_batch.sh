#!/usr/bin/env bash
# C15d: extract the ion-chamber summary for every run. Resumable.
#
#   ./ic_batch.sh [nparallel] [runlist]
#
# Per run: unpack ONLY the FRIB group -> a ~1 GB temp of raw traces, reduce it to the two numbers
# per event the analysis actually uses, then DELETE the temp.
#
#     icmax  = max of generic trace[0] over tb 1050-1250   (the ion chamber)
#     npulse = pulse count over tb 800-1500 above 200      (pile-up rejection)
#
# ★ NEVER KEEP THE TRACES. They are ~36 kB/event, i.e. ~100 GB over these 105 runs, and only those
# two numbers are ever read. The summary is ~130 kB per run -- ~14 MB for everything -- so the IC
# spectrum can be re-binned and the beam window re-chosen instantly, and the 100 GB never exists.
#
# The 1050-1250 window is inherited from a2091 and was CHECKED here before use: on run_0017 the
# mean trace peaks at tb 1140 and the per-event argmax sits at 1134-1141 across the interquartile
# range, so the window is centred on the pulse rather than clipping it.
#
# run_0047 has no /frib group at all and will fail; that is expected and it is skipped by name.

set -eo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"

NPAR="${1:-6}"
RUNLIST="${2:-$HERE/runs_d2.txt}"

IC_DIR="${C15D_IC:-/home/yassid/C15d_ic}"
TMP_DIR="$IC_DIR/tmp"
LOG_DIR="${C15D_LOGS:-/home/yassid/C15d_logs}"
MIN_FREE_GB="${MIN_FREE_GB:-40}"

set +u
# shellcheck disable=SC1091
source "$REPO/build/config.sh" >/dev/null 2>&1
set -u

[[ -n "${VMCWORKDIR:-}" ]] || { echo "ERROR: config.sh did not set VMCWORKDIR" >&2; exit 1; }
[[ -f "$RUNLIST" ]] || { echo "ERROR: run list not found: $RUNLIST" >&2; exit 1; }

mkdir -p "$IC_DIR" "$TMP_DIR" "$LOG_DIR"

mapfile -t RUNS < <(grep -vE '^\s*(#|$)' "$RUNLIST" | grep -v 'run_0047')
[[ ${#RUNS[@]} -gt 0 ]] || { echo "ERROR: empty run list -- refusing to no-op silently" >&2; exit 1; }

echo "=== C15d ic_batch ==="
echo "  runs     : ${#RUNS[@]} from $(basename "$RUNLIST")  (run_0047 skipped: no /frib)"
echo "  parallel : $NPAR"
echo "  out      : $IC_DIR"
echo

trap 'pkill -P $$ 2>/dev/null || true' EXIT INT TERM

do_ic() {
   local run="$1"
   local out="$IC_DIR/${run}_ic.root"
   local log="$LOG_DIR/${run}_ic.log"
   local tmp="$TMP_DIR/${run}_FRIB.root"

   if [[ -s "$out" ]]; then echo "[$run] ic exists, skipping"; return 0; fi

   local free_gb
   free_gb=$(df -BG --output=avail "$IC_DIR" | tail -1 | tr -d ' G')
   if (( free_gb < MIN_FREE_GB )); then echo "[$run] SKIP: only ${free_gb} GB free"; return 0; fi

   rm -f "$tmp"
   if ! root -b -q "$HERE/pid/unpackFRIB_C15d.C(\"$run\",\"$TMP_DIR/\")" >"$log" 2>&1 || [[ ! -s "$tmp" ]]; then
      rm -f "$tmp"
      # Distinguish "this run carries no ion-chamber data" from "the unpacker broke". Several runs
      # have a /frib group whose evt/ is EMPTY (0024, 0025) and one has no /frib at all (0047);
      # calling those a FAILURE sends the reader hunting for a bug that is not there.
      if ! grep -q "Unpacking FRIB for" "$log"; then
         echo "[$run] no IC data in this run (empty or absent /frib) -- skipped, not a failure"
      else
         echo "[$run] FRIB UNPACK FAILED -- see $log"
      fi
      return 0
   fi
   # Reduce, then delete the temp WHETHER OR NOT the reduction worked: a leftover 1 GB temp per
   # failed run fills the disk and stops every later run via the guard above.
   root -b -q "$HERE/pid/icsum_C15d.C(\"$run\",\"$TMP_DIR/\",\"$IC_DIR/\")" >>"$log" 2>&1 || true
   rm -f "$tmp"
   if [[ -s "$out" ]]; then
      echo "[$run] ic OK  ($(grep -oP 'peak=\K[0-9]+' "$log" | tail -1) ADC peak)"
   else
      echo "[$run] IC REDUCE FAILED -- see $log"
   fi
}
export -f do_ic
export HERE IC_DIR TMP_DIR LOG_DIR MIN_FREE_GB

printf '%s\n' "${RUNS[@]}" | xargs -P "$NPAR" -I{} bash -c 'do_ic "$@"' _ {}

echo
echo "=== done: $(ls -1 "$IC_DIR"/*_ic.root 2>/dev/null | wc -l) IC summaries, $(du -sh "$IC_DIR" 2>/dev/null | cut -f1) ==="
rmdir "$TMP_DIR" 2>/dev/null || true
