#!/usr/bin/env bash
# Stop the reconstruction cleanly after the run in flight, then validate the working point on
# whatever is already reconstructed.
#
#   ./validate_a2091_C15d.sh [nfit]
#
# WHY STOP EARLY. The reconstruction of 104 runs takes ~9 more hours, and every hour of it is
# spent on a drift working point that has NOT been checked against data. The a2091 par uses
# Spyral's dv = 1.1364 cm/us; the a2091 PROTON runs of the same experiment give a clean elastic
# ridge at 201 MeV with dv = 1.30. H2 and D2 at the same pressure have the same number density and
# hence the same E/N, so the two should agree closely -- they differ by 14 %. The (d,d) elastic
# ridge on these runs is the same beam and must also give ~201 MeV. If it does not, dv is wrong
# and every reconstructed run would have to be redone.
#
#   1  wait for the run in flight, then stop the reco loop
#   2  points file, no gain
#   3  per-run gain
#   4  points with gain + the PID plane
#   5  gain quality report          <- the user asked to see this before gating
#   6  GENFIT+CATIMA deuteron fits, UNGATED, on a subset
#   7  beam energy from the 15C(d,d) elastic ridge   <- the working-point test

set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
NFIT="${1:-8}"
LOG="${C15D_LOGS:-/home/yassid/C15d_logs}"
RECO="${C15D_RECO:-/home/yassid/C15d_reco}"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
[[ -n "${VMCWORKDIR:-}" ]] || { echo "ERROR: config.sh did not set VMCWORKDIR" >&2; exit 1; }
say(){ echo; echo "=== $(date '+%F %H:%M:%S')  $*"; }
cd "$HERE"

say "stage 1/7  waiting for the run in flight to finish, then stopping the reco"
n0=$(ls -1 "$RECO"/*_reco.root 2>/dev/null | wc -l)
# a run takes 150-900 s; give it 30 min, then stop regardless
for _ in $(seq 180); do
   n1=$(ls -1 "$RECO"/*_reco.root 2>/dev/null | wc -l)
   (( n1 > n0 )) && break
   pgrep -f "[u]npackReco_C15d" >/dev/null 2>&1 || break
   sleep 10
done
P1="production_""a2091"; P2="reco_chunked_""C15d.sh"
for p in $(ps -eo pid,cmd --no-headers | grep -F -e "$P1" -e "$P2" | grep -v grep | awk '{print $1}'); do
   kill -9 "$p" 2>/dev/null || true
done
sleep 2
for p in $(ps -eo pid,cmd --no-headers | grep '[r]oot.exe' | grep unpackReco | awk '{print $1}'); do
   kill -9 "$p" 2>/dev/null || true
done
# only NOW is .chunk safe to clear: nothing is writing to it. Clearing it while a run is live
# destroys finished output that has not been mv-ed yet -- that already cost two full runs once.
sleep 3
rm -rf "$RECO"/.chunk/* 2>/dev/null || true
NR=$(ls -1 "$RECO"/*_reco.root 2>/dev/null | wc -l)
echo "  reconstruction stopped with $NR runs done"
[[ "$NR" -ge 10 ]] || { echo "  too few runs to validate -- stopping"; exit 1; }

say "stage 2/7  points file, no gain"
root -b -q 'pid/make_points_C15d.C("/home/yassid/C15d_reco/","pid/points_C15d.root","")' \
   >>"$LOG/val_points.log" 2>&1 || true

say "stage 3/7  measure the per-run gain"
root -b -q 'measure_gain_C15d.C(0.30,0.40,"/home/yassid/C15d_reco/","gainmatch_C15d.csv")' \
   >>"$LOG/val_gain.log" 2>&1 || true
grep -E "estimator|reference|factor range|interpolated" "$LOG/val_gain.log" | tail -4 || true

say "stage 4/7  points with gain, and the PID plane"
root -b -q 'pid/make_points_C15d.C("/home/yassid/C15d_reco/","pid/points_C15d.root","gainmatch_C15d.csv")' \
   >>"$LOG/val_points.log" 2>&1 || true
root -b -q 'mkpid_C15d.C("/home/yassid/C15d_reco/","plots/",340,0,85,300,0,2.5,"",0,-1,"gainmatch_C15d.csv")' \
   >>"$LOG/val_plane.log" 2>&1 || true
grep -E "runs|tracks|selected|wrote" "$LOG/val_plane.log" | tail -4 || true

say "stage 5/7  GAIN QUALITY REPORT"
root -b -q 'gain_report_C15d.C()' 2>&1 | sed 's/\x1b\[[0-9;]*m//g' | sed -n '/GAIN QUALITY/,$p' | tee "$LOG/val_gainreport.log"

say "stage 6/7  GENFIT+CATIMA deuteron fits, ungated, $NFIT runs"
ls -S "$RECO"/*_reco.root 2>/dev/null | head -"$NFIT" | xargs -r -n1 basename | sed 's/_reco.root//' > "$HERE/.val_runs.txt"
./fit_batch.sh 8 d "$HERE/.val_runs.txt" -1 "" kTRUE >>"$LOG/val_fit.log" 2>&1 || true
echo "  kin files: $(ls -1 /home/yassid/C15d_fit/*_kin_d.root 2>/dev/null | wc -l)"

say "stage 7/7  beam energy from the 15C(d,d) elastic ridge"
root -b -q 'dd/kin_dd_C15d.C("/home/yassid/C15d_fit/","","dd/plots/")' >>"$LOG/val_eb.log" 2>&1 || true
# Scan the KE floor rather than trusting one value: it is CHANNEL-SPECIFIC and the default here
# was tuned on the a1975 16C data. Too high a floor cuts away the very ridge it is meant to find
# (that produced "no ridge in the seed slice" on the C15p protons); too low and the finder locks
# onto the low-energy bulk (that gave 20 MeV instead of 188 on a1975).
for f in 3 8 12 16 20 24; do
   printf "  keMinUse %-3s : " "$f"
   root -b -q "dd/ebeam_dd_C15d.C(\"dd/plots/dd_kin_C15d.root\",\"dd/plots/\",42,78,2.0,50.0,$f)" 2>&1 \
     | sed 's/\x1b\[[0-9;]*m//g' | grep -oE "Ebeam = [0-9]+ MeV, window-to-window spread [0-9-]+ MeV|ERROR.*" | head -1
done | tee -a "$LOG/val_eb.log"
root -b -q 'dd/ebeam_dd_C15d.C("dd/plots/dd_kin_C15d.root","dd/plots/")' 2>&1 \
   | sed 's/\x1b\[[0-9;]*m//g' | grep -E "locus-count" | tee -a "$LOG/val_eb.log"
echo
echo "  THE TEST: the a2091 PROTON runs of this same experiment give 201 MeV (13.4 MeV/u)."
echo "  Same beam, so the (d,d) ridge must agree. A large disagreement means dv 1.1364 is wrong."

say "VALIDATION DONE"
touch "$HERE/.validation_DONE"
