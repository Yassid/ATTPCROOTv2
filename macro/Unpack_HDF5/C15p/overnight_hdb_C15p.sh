#!/usr/bin/env bash
# C15p overnight chain, MultiFit + HDBSCAN rebuild. Runs everything that does NOT need a human,
# and STOPS at the PID plane so the gates can be drawn on it in the morning.
#
#   ./overnight_hdb_C15p.sh [nparallel]
#
#   1  full reco from raw HDF5   AtPSAMultiFit + AtDirDeDxCleaner + HDBSCAN, doSC=false,
#                                AtPIDTask in-chain  -> AtPIDEvent is already persisted
#   2  flatten AtPIDEvent to per-run _pid.root ntuples
#   3  points file, NO gain
#   4  per-run gain from that plane
#   5  points again WITH gain, and the PID plane PNG   <- draw the gates on this
#
# There is NO pidpass stage: stage 1 writes AtPIDEvent itself, which the old chain could not do
# because it was reusing a reco that predated the branch.
# There is NO fitting stage: fitting before the gates exist is what cost 43 GB last time.
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
NPAR="${1:-24}"
RECO="${C15P_RECO_HDB:-/home/yassid/Data/a2091_C15p_reco_hdb}"
LOG="${C15P_LOGS_HDB:-/home/yassid/C15p_hdb_logs}"
RUNLIST="$HERE/runs_pp_hdb.txt"
mkdir -p "$LOG" "$HERE/plots_hdb" "$HERE/pid"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
[[ -n "${VMCWORKDIR:-}" ]] || { echo "ERROR: config.sh did not set VMCWORKDIR" >&2; exit 1; }
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
cd "$HERE"
say(){ echo; echo "=== $(date '+%F %H:%M:%S')  $*"; }
NRUNS=$(grep -vE '^\s*(#|$)' "$RUNLIST" | tr -s ' \n' '\n' | grep -c '^run_')

say "stage 1/5  reco from raw: MultiFit + HDBSCAN, space charge OFF  ($NRUNS runs, $NPAR parallel)"
./reco_batch_C15p.sh "$NPAR" "$RUNLIST" >>"$LOG/reco.log" 2>&1 || true
NR=$(ls -1 "$RECO"/*_reco.root 2>/dev/null | wc -l)
echo "  reco files: $NR / $NRUNS"
[[ "$NR" -ge 10 ]] || { echo "  too few to continue -- stopping"; exit 1; }

say "stage 2/5  flatten AtPIDEvent -> per-run _pid.root"
flat(){ local r="$1"
  [[ -s "$RECO/${r}_pid.root" ]] && return 0
  [[ -s "$RECO/${r}_reco.root" ]] || return 0
  root -b -q "$HERE/pidntuple_C15p.C(\"$r\",\"$RECO/\",\"$RECO/\")" >>"$LOG/${r}_pidntuple.log" 2>&1 || true
}
export -f flat; export RECO HERE LOG
grep -vE '^\s*(#|$)' "$RUNLIST" | tr -s ' \n' '\n' | grep '^run_' \
  | xargs -P "$NPAR" -I{} bash -c 'flat "$@"' _ {}
echo "  pid ntuples: $(ls -1 "$RECO"/*_pid.root 2>/dev/null | wc -l) / $NRUNS"

say "stage 3/5  points file, no gain"
root -b -q "pid/make_points_C15p.C(\"$RECO/\",\"pid/points_hdb_C15p.root\",\"\")" \
   >>"$LOG/points_nogain.log" 2>&1 || true

say "stage 4/5  measure the per-run gain"
root -b -q "measure_gain_C15p.C(0.30,0.40,\"$RECO/\",\"gainmatch_hdb_C15p.csv\")" \
   >>"$LOG/gain.log" 2>&1 || true
grep -E "estimator|reference|factor range" "$LOG/gain.log" | tail -3 || true

say "stage 5/5  points WITH gain, and the PID plane for the gates"
root -b -q "pid/make_points_C15p.C(\"$RECO/\",\"pid/points_hdb_C15p.root\",\"gainmatch_hdb_C15p.csv\")" \
   >>"$LOG/points_gain.log" 2>&1 || true
root -b -q "mkpid_C15p.C(\"$RECO/\",\"plots_hdb/\",340,0,85,300,0,2.5,\"\",0,-1,\"gainmatch_hdb_C15p.csv\")" \
   >>"$LOG/plane.log" 2>&1 || true
grep -E "tracks|selected|wrote" "$LOG/plane.log" | tail -3 || true

say "DONE -- draw the gates on plots_hdb/pid_C15p.png  (pid/open_gate_draw.sh)"
touch "$HERE/.overnight_hdb_DONE"
