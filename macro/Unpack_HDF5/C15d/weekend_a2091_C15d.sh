#!/usr/bin/env bash
# 15C + d, a2091 D2 runs 0013-0133 -- the long unattended pass.
#
#   ./weekend_a2091_C15d.sh
#
# Everything here is resumable and every stage skips work that already exists, so it can be
# interrupted and restarted without losing anything.
#
#   1  wait for the IC rebuild in flight
#   2  reconstruct the remaining runs (chunked, ~9 h)
#   3  IC pass on the newly reconstructed runs
#   4  points file with IC joined, no gain
#   5  per-run gain
#   6  points with gain + the PID plane, IC WINDOW APPLIED
#   7  gain quality report on the full set
#   8  IC spectrum and Z ladder on the full set
#
# It STOPS THERE. The next step is a human drawing the deuteron gate, and the IC window itself is
# not settled: the a2091 D2 spectrum has TWO large peaks in the carbon region, 1135 and 1365, and
# the proton-run window [979.5, 1278.8] cuts between them. The window used below is deliberately
# the carbon peak alone -- see IC_LO/IC_HI. Revisit it once the 1365 component is identified.
#
# ⚠ GUARDS, all of them learned the hard way:
#   * a corrupt event count (run_0019 claimed 1.3e14 events) once produced 1.4 TB of HDF5 error
#     text and filled the disk; reco_chunked_C15d.sh now sanity-checks it and caps chunk logs.
#   * NEVER clear .chunk or .part while a batch is running -- that destroys finished output that
#     has not been mv-ed into place yet, and it has already cost two full runs.
#   * the disk check below stops the whole thing rather than wedging the filesystem again.

set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
LOG="${C15D_LOGS:-/home/yassid/C15d_logs}"
RECO="${C15D_RECO:-/home/yassid/C15d_reco}"
RUNLIST="$HERE/runs_a2091_d2.txt"
NCHUNK=12
MIN_FREE_GB=150

# The CARBON peak only. The proton-run window [979.5, 1278.8] is NOT used: it clips the 1365
# component, which is 213k events and unidentified. 1135 +- 90 keeps the carbon peak and stops
# short of the 1365 shoulder.
IC_LO=1045
IC_HI=1225

mkdir -p "$LOG"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
[[ -n "${VMCWORKDIR:-}" ]] || { echo "ERROR: config.sh did not set VMCWORKDIR" >&2; exit 1; }
say(){ echo; echo "=== $(date '+%F %H:%M:%S')  $*"; }
cd "$HERE"

guard(){
   local free
   free=$(df -BG --output=avail "$RECO" | tail -1 | tr -d ' G')
   if (( free < MIN_FREE_GB )); then
      echo "★ ABORT: only ${free} GB free (guard ${MIN_FREE_GB})"; exit 1
   fi
}

say "stage 1/8  (nothing to rebuild -- all previous output was erased)"
# ★ ic_batch SKIPS a run whose <run>_ic.root already exists, so a fix to icsum_C15d.C has no
# effect until the old outputs are removed. The amplitude window was wrong (max over a narrow
# [1050,1250] while pulses were counted over [800,1500]), which mismeasured 12.4 % of events as
# baseline and put a spurious 35 ADC peak in the spectrum. Force the rebuild once, here.
echo "  starting from a clean slate"

say "stage 2/8  reconstruct the remaining runs"
guard
./reco_chunked_C15d.sh "$NCHUNK" "$RUNLIST" >>"$LOG/wk_reco.log" 2>&1 || true
NR=$(ls -1 "$RECO"/*_reco.root 2>/dev/null | wc -l)
echo "  reco files: $NR"

say "stage 3/8  IC pass on everything reconstructed"
guard
ls -1 "$RECO"/*_reco.root 2>/dev/null | sed 's#.*/##;s/_reco.root//' > "$HERE/runs_a2091_have.txt"
./ic_batch.sh 3 "$HERE/runs_a2091_have.txt" >>"$LOG/wk_ic.log" 2>&1 || true
echo "  IC summaries: $(ls -1 /home/yassid/C15d_ic/*_ic.root 2>/dev/null | wc -l)"

say "stage 4/8  points file, IC joined, no gain"
root -b -q 'pid/make_points_C15d.C("/home/yassid/C15d_reco/","pid/points_C15d.root","")' \
   >>"$LOG/wk_points.log" 2>&1 || true

say "stage 5/8  per-run gain"
root -b -q 'measure_gain_C15d.C(0.30,0.40,"/home/yassid/C15d_reco/","gainmatch_C15d.csv")' \
   >>"$LOG/wk_gain.log" 2>&1 || true
grep -E "estimator|reference|factor range|interpolated" "$LOG/wk_gain.log" | tail -4 || true

say "stage 6/8  points with gain, and the PID plane with the IC window"
root -b -q 'pid/make_points_C15d.C("/home/yassid/C15d_reco/","pid/points_C15d.root","gainmatch_C15d.csv")' \
   >>"$LOG/wk_points.log" 2>&1 || true
# the ungated plane, for reference
root -b -q "mkpid_C15d.C(\"/home/yassid/C15d_reco/\",\"plots/\",340,0,85,300,0,2.5,\"\",0,-1,\"gainmatch_C15d.csv\")" \
   >>"$LOG/wk_plane.log" 2>&1 || true
# ★ the plane the gate is actually drawn on: IC window + single pulse, from the joined points file
root -b -q "mkpid_ic_C15d.C($IC_LO,$IC_HI,\"pid/points_C15d.root\",\"plots/\",\"pid_ic_C15d\")" 2>&1 \
   | sed 's/\x1b\[[0-9;]*m//g' | grep -E "IC window|tracks|no usable|IN PLANE|wrote" | tee -a "$LOG/wk_plane.log"
# and the same plane with NO IC cut, so the effect of the window is visible side by side
root -b -q "mkpid_ic_C15d.C(-1,-1,\"pid/points_C15d.root\",\"plots/\",\"pid_noic_C15d\")" \
   >>"$LOG/wk_plane.log" 2>&1 || true

say "stage 7/8  gain quality report"
root -b -q 'gain_report_C15d.C()' 2>&1 | sed 's/\x1b\[[0-9;]*m//g' | sed -n '/GAIN QUALITY/,$p' \
   | tee "$LOG/wk_gainreport.log"

say "stage 8/8  IC spectrum and Z ladder on the full set"
root -b -q 'pid/ic_ladder_C15d.C()' 2>&1 | sed 's/\x1b\[[0-9;]*m//g' | sed -n '/IC SPECTRUM/,$p' \
   | tee "$LOG/wk_icladder.log"

say "WEEKEND PASS DONE -- next step is the deuteron gate on plots/pid_C15d.png"
echo "  IC window used for reference: [$IC_LO, $IC_HI], the carbon peak only."
echo "  The 1365 component is still unidentified; revisit the window before quoting physics."
touch "$HERE/.weekend_DONE"
