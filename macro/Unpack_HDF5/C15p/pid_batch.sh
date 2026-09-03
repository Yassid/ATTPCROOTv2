#!/usr/bin/env bash
# C15p: add an AtPIDEvent branch to the existing a2091 reconstruction.
#
#   ./pid_batch.sh [nparallel] [runlist]
#
# The a2091 reco was built by the older ATTPCROOTv2 tree and carries AtEventH + AtPatternEvent
# but NO AtPIDEvent. This runs AtPIDTask over it and writes a small companion tree to
# $C15P_RECO. Resumable: a run whose output already exists is skipped.
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
NPAR="${1:-16}"
RUNLIST="${2:-$HERE/runs_pp.txt}"
OUT="${C15P_RECO:-/home/yassid/C15p_reco}"
LOG="${C15P_LOGS:-/home/yassid/C15p_logs}"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
[[ -n "${VMCWORKDIR:-}" ]] || { echo "ERROR: config.sh did not set VMCWORKDIR" >&2; exit 1; }
mkdir -p "$OUT" "$LOG"
cd "$HERE"
mapfile -t RUNS < <(grep -vE '^\s*(#|$)' "$RUNLIST" | tr -s ' \n' '\n' | grep '^run_')
echo "=== C15p pid_batch : ${#RUNS[@]} runs, $NPAR parallel -> $OUT ==="
one(){
  local run="$1"
  if [[ -s "$OUT/${run}_pid.root" ]]; then echo "[$run] pid exists, skipping"; return 0; fi
  if [[ ! -s "/home/yassid/a2091_C15_reco/${run}_reco.root" ]]; then echo "[$run] NO reco, skipping"; return 0; fi
  # TWO steps, and both are needed: pidpass writes the AtPIDEvent branch, pidntuple flattens it
  # to <run>_pid.root, which is what make_points_C15p and mkpid_C15p actually read. Stopping
  # after pidpass leaves the plane stage with nothing to chain and it fails silently.
  if ! root -b -q "$HERE/pidpass_C15p.C(\"$run\",-1)" >"$LOG/${run}_pid.log" 2>&1 \
     || [[ ! -s "$OUT/${run}_reco.root" ]]; then
     echo "[$run] pidpass FAILED (see $LOG/${run}_pid.log)"; return 0
  fi
  if root -b -q "$HERE/pidntuple_C15p.C(\"$run\",\"$OUT/\")" >>"$LOG/${run}_pid.log" 2>&1 \
     && [[ -s "$OUT/${run}_pid.root" ]]; then
     echo "[$run] pid OK ($(du -h "$OUT/${run}_pid.root" | cut -f1) ntuple)"
  else
     echo "[$run] pidntuple FAILED (see $LOG/${run}_pid.log)"
  fi
}
export -f one; export OUT LOG HERE
printf '%s\n' "${RUNS[@]}" | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "=== done: $(ls -1 "$OUT"/*_reco.root 2>/dev/null | wc -l) pid files ==="
