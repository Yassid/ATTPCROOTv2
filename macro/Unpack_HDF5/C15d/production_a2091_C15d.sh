#!/usr/bin/env bash
# 15C + d production on the a2091 DEUTERIUM runs 0013-0133.
#
#   ./production_a2091_C15d.sh [nparallel]
#
# Stops after the PID plane ON PURPOSE. The next step is a human drawing the deuteron gate on
# plots/pid_C15d.png, and the user asked specifically to see the GAIN MATCHING first, since a
# gate drawn on an unmatched plane selects the wrong tracks run by run.
#
#   1  reco, all 104 runs                        (resumable; skips runs already done)
#   2  points file, NO gain                      -> the plane the gain is measured on
#   3  per-run gain factors                      -> gainmatch_C15d.csv  + diagnostics
#   4  points again WITH gain, and the PID plane -> the gate gets drawn on this
#   5  gain-quality report                       -> is the match good enough to gate on?
#
# ⚠ THIS IS a2091, NOT a1975. a1975 runs 17-103 are 16C(d,p)17C; the earlier C15d analysis read
# them by mistake and its 1.8 MeV Ex peak was the 16C 2+ at 1.766 MeV. See
# parameters/ATTPC.C15d_a2091_D2.par for the full account.
#
# ⚠ The unpacker prints "Number of events from metaData (N) does not match the number of entries
# in HDF5 file (2N)" on every a2091 file. It is SPURIOUS: /get holds two datasets per event
# (evtN_data and evtN_header), so the child count is twice the event count. The metadata number
# is the right one and is what the unpacker uses.

set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
NPAR="${1:-12}"   # chunks per run, not concurrent runs
LOG="${C15D_LOGS:-/home/yassid/C15d_logs}"
# Second argument so a production can be started alongside a batch that is still
# finishing: reco_batch stages into a PER-RUN .part directory, so two batches on
# DIFFERENT runs are safe, and stages 2-5 read the whole reco dir regardless.
RUNLIST="${2:-$HERE/runs_a2091_d2.txt}"
mkdir -p "$LOG"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
[[ -n "${VMCWORKDIR:-}" ]] || { echo "ERROR: config.sh did not set VMCWORKDIR" >&2; exit 1; }
say(){ echo; echo "=== $(date '+%F %H:%M:%S')  $*"; }
cd "$HERE"

say "stage 1/5  reco, $(grep -c '^run_' "$RUNLIST") runs from $(basename "$RUNLIST")"
# CHUNKED, not reco_batch: the raw data is on one spinning USB disk, so reconstructing several
# runs at once collapses it to seeks (12 runs -> 4.0 MB/s, 16 kB per I/O, 16 % CPU each). Reading
# ONE run into page cache and splitting it across cores instead gives 96 % CPU on 11.6 cores with
# the disk at 0.0 MB/s. Measured on run_0110, 24929 events: 3278 s -> 190 s, a 17x speedup.
# NPAR here is CHUNKS PER RUN, capped by memory (~1.1 GB per chunk + a 17-27 GB warm file in
# 62 GB of RAM), not by core count.
./reco_chunked_C15d.sh "$NPAR" "$RUNLIST" >>"$LOG/prod_reco.log" 2>&1 || true
NR=$(ls -1 /home/yassid/C15d_reco/*_reco.root 2>/dev/null | wc -l)
echo "  reco files: $NR"
[[ "$NR" -ge 20 ]] || { echo "  too few recos to continue -- stopping"; exit 1; }

say "stage 2/5  points file, no gain"
root -b -q 'pid/make_points_C15d.C("/home/yassid/C15d_reco/","pid/points_C15d.root","")' \
   >>"$LOG/prod_points.log" 2>&1 || true

say "stage 3/5  measure the per-run gain"
root -b -q 'measure_gain_C15d.C(0.30,0.40,"/home/yassid/C15d_reco/","gainmatch_C15d.csv")' \
   >>"$LOG/prod_gain.log" 2>&1 || true
grep -E "estimator|reference|factor range|interpolated" "$LOG/prod_gain.log" | tail -5 || true

say "stage 4/5  points WITH gain, and the PID plane"
root -b -q 'pid/make_points_C15d.C("/home/yassid/C15d_reco/","pid/points_C15d.root","gainmatch_C15d.csv")' \
   >>"$LOG/prod_points.log" 2>&1 || true
root -b -q 'mkpid_C15d.C("/home/yassid/C15d_reco/","plots/",340,0,85,300,0,2.5,"",0,-1,"gainmatch_C15d.csv")' \
   >>"$LOG/prod_plane.log" 2>&1 || true
grep -E "runs|tracks|selected|wrote" "$LOG/prod_plane.log" | tail -4 || true

say "stage 5/5  gain-quality report"
root -b -q 'gain_report_C15d.C()' >>"$LOG/prod_gainreport.log" 2>&1 || true
sed 's/\x1b\[[0-9;]*m//g' "$LOG/prod_gainreport.log" | sed -n '/GAIN QUALITY/,$p' || true

say "DONE -- next step is the DEUTERON GATE on plots/pid_C15d.png"
touch "$HERE/.production_DONE"
