#!/usr/bin/env bash
# C15p: regenerate the ion-chamber summaries WITH the true event number (evtid).
#
#   ./ic_batch_C15p.sh [nparallel] [runlist]
#
# The shipped IC files carry only `entry`, the position in the FRIB tree. /get and /frib are
# separate DAQs in one merged file: they do not record the same set, and a single empty event
# offsets every later entry (run_0140: 32,295 entries, evtid 0..32,293, one empty). A positional
# join therefore mismatches everything after the first gap, which is why runs 0140 and 0151 were
# refused outright and why only 7,059 of 18,264 fitted tracks could be matched to a beam value.
#
# icsum_C15p.C now writes `evtid` from AtRawEvent::GetEventID, and also takes icmax over the SAME
# range the pulses are counted over -- the narrow [1050,1250] window returned the BASELINE for any
# event whose pulse fell outside it while npulse still said 1, measured at 12.9 % on these runs.
#
# Two stages per run: unpackFRIB (raw HDF5 -> <run>_FRIB.root) then icsum -> <run>_ic.root.
# The FRIB intermediate is ~40 kB/event and is DELETED once summarised, so the staging cost is
# bounded by the parallelism rather than by the run set.
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
NPAR="${1:-4}"
RUNLIST="${2:-$HERE/runs_pp_hdb.txt}"
RAW="/media/yassid/Seagate Hub/ATTPC/Data/a2091/"
TMP="${C15P_IC_TMP:-/home/yassid/Data/a2091_C15p_ic_tmp}"
OUT="${C15P_IC:-/home/yassid/Data/a2091_C15p_ic_evtid}"
LOG="${C15P_LOGS:-/home/yassid/C15p_hdb_logs}/ic"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
[[ -d "$RAW" ]] || { echo "ERROR: raw not mounted: $RAW" >&2; exit 1; }
mkdir -p "$TMP" "$OUT" "$LOG"; cd "$HERE"
mapfile -t RUNS < <(grep -vE '^\s*(#|$)' "$RUNLIST" | tr -s ' \n' '\n' | grep '^run_')
echo "=== C15p IC regeneration: ${#RUNS[@]} runs, $NPAR parallel -> $OUT ==="

one(){
  local r="$1"
  if [[ -s "$OUT/${r}_ic.root" ]]; then echo "[$(date +%H:%M:%S)] $r  ic exists, skipping"; return 0; fi
  [[ -f "$RAW/${r}.h5" ]] || { echo "[$(date +%H:%M:%S)] $r  NO RAW"; return 0; }
  local t0=$SECONDS
  if [[ ! -s "$TMP/${r}_FRIB.root" ]]; then
     root -b -q "$HERE/pid/unpackFRIB_C15p.C(\"$r\",\"$TMP/\",\"$RAW\")" >"$LOG/${r}_frib.log" 2>&1 || true
  fi
  [[ -s "$TMP/${r}_FRIB.root" ]] || { echo "[$(date +%H:%M:%S)] $r  FRIB FAILED (see $LOG/${r}_frib.log)"; return 0; }
  root -b -q "$HERE/pid/icsum_C15p.C(\"$r\",\"$TMP/\",\"$OUT/\")" >"$LOG/${r}_icsum.log" 2>&1 || true
  if [[ -s "$OUT/${r}_ic.root" ]]; then
     echo "[$(date +%H:%M:%S)] $r  OK  $(grep -oE 'entries=[0-9]+ +withIC=[0-9]+' "$LOG/${r}_icsum.log" | head -1)  in $(( (SECONDS-t0)/60 )) min"
     rm -f "$TMP/${r}_FRIB.root"          # bounded staging: summarised, so the intermediate goes
  else
     echo "[$(date +%H:%M:%S)] $r  ICSUM FAILED (see $LOG/${r}_icsum.log)"
  fi
}
export -f one; export RAW TMP OUT LOG HERE
printf '%s\n' "${RUNS[@]}" | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "=== done: $(ls -1 "$OUT"/*_ic.root 2>/dev/null | wc -l) / ${#RUNS[@]} ic files ==="
