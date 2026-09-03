#!/usr/bin/env bash
# C15p overnight chain: 15C(p,p') on the a2091 proton-target runs, GENFIT + CATIMA.
#
#   ./overnight_C15p.sh [nparallel]
#
# There is NO reconstruction stage: the 41 a2091 runs were already reconstructed by the older
# ATTPCROOTv2 tree (86 GB) and those files were verified readable by THIS build before the port
# (run_0138: 333 tracks / 95381 hits over 300 events), so the expensive stage is skipped entirely.
# run_0140 is excluded -- its reco is a 288-byte truncated file.
#
#   1  PID pass over the existing reco            -> AtPIDEvent
#   2  points file (IC joined), no gain yet
#   3  per-run gain from that plane
#   4  points again with gain + the PID plane PNG  <- the proton gate gets drawn on this
#   5  GENFIT+CATIMA proton fits, UNGATED          <- the beam-energy check
#   6  beam energy from the 15C(p,p) elastic ridge
#
# Stage 5 is ungated because the proton gate needs a human. The elastic ridge does not need a PID
# gate, so the beam energy can still be checked overnight -- and it needs checking: an ungated
# fit of run_0138 gave 196 MeV against the 157 MeV the old a2091 analysis quotes.
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
NPAR="${1:-16}"
LOG="${C15P_LOGS:-/home/yassid/C15p_logs}"
mkdir -p "$LOG"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
[[ -n "${VMCWORKDIR:-}" ]] || { echo "ERROR: config.sh did not set VMCWORKDIR" >&2; exit 1; }
say(){ echo; echo "=== $(date '+%H:%M:%S')  $*"; }
cd "$HERE"

say "stage 1/6  PID pass over the existing a2091 reco"
./pid_batch.sh "$NPAR" >>"$LOG/overnight_pid.log" 2>&1 || true
NP=$(ls -1 /home/yassid/C15p_reco/*_reco.root 2>/dev/null | wc -l)
echo "  pid files: $NP / 40"
[[ "$NP" -ge 10 ]] || { echo "  too few to continue -- stopping"; exit 1; }

say "stage 2/6  points file, no gain"
root -b -q 'pid/make_points_C15p.C("/home/yassid/C15p_reco/","pid/points_C15p.root","")' \
   >>"$LOG/overnight_points.log" 2>&1 || true

say "stage 3/6  measure the gain"
root -b -q 'measure_gain_C15p.C(0.30,0.40,"/home/yassid/C15p_reco/","gainmatch_C15p.csv")' \
   >>"$LOG/overnight_gain.log" 2>&1 || true
grep -E "estimator|reference|factor range" "$LOG/overnight_gain.log" | tail -3 || true

say "stage 4/6  points with gain, and the PID plane"
root -b -q 'pid/make_points_C15p.C("/home/yassid/C15p_reco/","pid/points_C15p.root","gainmatch_C15p.csv")' \
   >>"$LOG/overnight_points.log" 2>&1 || true
root -b -q 'mkpid_C15p.C("/home/yassid/C15p_reco/","plots/",340,0,85,300,0,2.5,"",0,-1,"gainmatch_C15p.csv")' \
   >>"$LOG/overnight_plane.log" 2>&1 || true
grep -E "tracks|selected|wrote" "$LOG/overnight_plane.log" | tail -3 || true

say "stage 5/6  GENFIT+CATIMA proton fits, ungated, all 40 runs"
./fit_batch.sh "$NPAR" p "$HERE/runs_pp.txt" -1 "" kTRUE >>"$LOG/overnight_fit.log" 2>&1 || true
echo "  kin files: $(ls -1 /home/yassid/C15p_fit/*_kin_p.root 2>/dev/null | wc -l)"

say "stage 6/6  beam energy from the elastic ridge"
root -b -q 'pp/kin_pp_C15p.C("/home/yassid/C15p_fit/","","pp/plots/")' >>"$LOG/overnight_eb.log" 2>&1 || true
root -b -q 'pp/ebeam_pp_C15p.C("pp/plots/pp_kin_C15p.root","pp/plots/")' >>"$LOG/overnight_eb.log" 2>&1 || true
grep -E "locus-count|Ebeam =" "$LOG/overnight_eb.log" | tail -4 || true

say "DONE"
touch "$HERE/.overnight_DONE"
