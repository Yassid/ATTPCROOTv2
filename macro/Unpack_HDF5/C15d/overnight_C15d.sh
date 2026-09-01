#!/usr/bin/env bash
# C15d overnight chain at working point A (dv 0.556, 3.125 MHz, TBentrance 582).
#
#   ./overnight_C15d.sh [nparallel]
#
# Stages, each resumable and each skipped if its output already exists:
#   1  reco + PID cache, all 105 runs                  ~10-12 h at 8-parallel
#   2  points file (IC joined, no gain yet)
#   3  measure the per-run gain from this plane
#   4  points file again, gain applied, and the PID plane PNG   <- gates get drawn on this
#   5  GENFIT+CATIMA deuteron fits, UNGATED, on a subset        <- the beam-energy check
#   6  beam energy from the 15C(d,d) elastic ridge
#
# Stage 5 is ungated because the gates must be REDRAWN on the new plane and that needs a human.
# The elastic ridge does not need a PID gate, so the beam energy can still be checked overnight.
#
# The fitter settings are the ones the 14C(d,p) simulation established: backwardSeedFix ON,
# rangeConstraint OFF (it reaches 4 % of backward tracks and destroys 9 % of forward ones),
# material effects ON with the CATIMA backends, matFallback OFF.

set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
NPAR="${1:-8}"
LOG=/home/yassid/C15d_logs
mkdir -p "$LOG"

set +u
# shellcheck disable=SC1091
source "$REPO/build/config.sh" >/dev/null 2>&1
set -u
[[ -n "${VMCWORKDIR:-}" ]] || { echo "ERROR: config.sh did not set VMCWORKDIR" >&2; exit 1; }

say() { echo; echo "=== $(date '+%H:%M:%S')  $*"; }
cd "$HERE"

say "stage 1/6  reco + PID cache, 105 runs"
./reco_batch.sh "$NPAR" >>"$LOG/overnight_reco.log" 2>&1 || true
NRECO=$(ls -1 /home/yassid/C15d_reco/*_reco.root 2>/dev/null | wc -l)
echo "  recos: $NRECO / 105"
[[ "$NRECO" -ge 20 ]] || { echo "  too few recos to continue -- stopping"; exit 1; }

say "stage 2/6  points file, no gain"
root -b -q 'pid/make_points_C15d.C("/home/yassid/C15d_reco/","pid/points_C15d.root","")' \
   >>"$LOG/overnight_points.log" 2>&1 || true

say "stage 3/6  measure the gain"
root -b -q 'measure_gain_C15d.C(0.30,0.40,"/home/yassid/C15d_reco/","gainmatch_C15d.csv")' \
   >>"$LOG/overnight_gain.log" 2>&1 || true
grep -E "estimator|reference|measured |factor range" "$LOG/overnight_gain.log" | tail -4 || true

say "stage 4/6  points file with gain, and the PID plane"
root -b -q 'pid/make_points_C15d.C("/home/yassid/C15d_reco/","pid/points_C15d.root","gainmatch_C15d.csv")' \
   >>"$LOG/overnight_points.log" 2>&1 || true
root -b -q 'mkpid_C15d.C("/home/yassid/C15d_reco/","plots/",340,0,85,300,0,2.5,"",0,-1,"gainmatch_C15d.csv")' \
   >>"$LOG/overnight_plane.log" 2>&1 || true
grep -E "tracks|selected|wrote" "$LOG/overnight_plane.log" | tail -3 || true

say "stage 5/6  GENFIT+CATIMA deuteron fits, ungated, subset"
# D2 ONLY. Sorting by size and taking the top 12 selects the HYDROGEN runs (>=106), which are
# the largest files on disk -- 7 of the 12 came out H2 the first time this ran, and fitting a
# D2/H2 mixture as one channel put the (d,d) elastic ridge at 54 MeV instead of 190. Filter on
# the run number FIRST, then take the largest of what remains.
ls -S /home/yassid/C15d_reco/*_reco.root 2>/dev/null \
   | xargs -r -n1 basename | sed 's/_reco.root//' \
   | awk -F_ '$2+0 <= 103' | head -12 > "$HERE/.overnight_runs.txt"
while read -r r; do
   [[ -n "$r" ]] || continue
   [[ -s "/home/yassid/C15d_fit/${r}_kin_d.root" ]] && continue
   root -b -q "fitGenfit_C15d.C(\"$r\", 12000, \"/home/yassid/C15d_reco/\", \"\", \"/home/yassid/C15d_fit/\", -2.85, 2, 5, \"\", 4.0, 5.0, 178.0, kTRUE, kTRUE, 1000010020, 2.01410178, 1, \"d\")" \
      >>"$LOG/overnight_fit.log" 2>&1 || true
   root -b -q "dump_kine_C15d.C(\"$r\",\"\",\"/home/yassid/C15d_fit/\",\"\",\"d\")" \
      >>"$LOG/overnight_fit.log" 2>&1 || true
done < "$HERE/.overnight_runs.txt"
echo "  kin files: $(ls -1 /home/yassid/C15d_fit/*_kin_d.root 2>/dev/null | wc -l)"

say "stage 6/6  beam energy from the elastic ridge"
root -b -q 'dd/kin_dd_C15d.C("/home/yassid/C15d_fit/","","dd/plots/")' >>"$LOG/overnight_eb.log" 2>&1 || true
root -b -q 'dd/ebeam_dd_C15d.C()' >>"$LOG/overnight_eb.log" 2>&1 || true
grep -E "locus-count|Ebeam =" "$LOG/overnight_eb.log" | tail -4 || true

say "DONE"
touch "$HERE/.overnight_DONE"
