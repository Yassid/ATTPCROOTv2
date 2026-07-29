#!/bin/bash
# ============================================================================
# Overnight chain for a2091 15C+p, run after the Phase B reco + first fit pass
# were already launched. Waits for each stage rather than assuming timings.
#
#   S1  wait for the FRIBDAQ reco (overnight_C15.sh Phase B) to finish
#   S2  wait for the first fitall pass (28 standard runs) to finish
#   S3  fitall again -> picks up the 11 newly-reco'd FRIBDAQ runs
#   S4  physics: Ex spectra per fitter + Ebeam scan + fitter comparison
#
# Everything in S3/S4 is idempotent, so re-running the chain is safe.
# ============================================================================
REPO=/home/yassid/fair_install/ATTPCROOTv2
HERE=$REPO/macro/Unpack_HDF5/a2091/UKF
RECO=/home/yassid/a2091_C15_reco
FIT=/home/yassid/a2091_C15_fit
LOG=$RECO/logs
CHAIN=$RECO/chain.log
mkdir -p "$LOG"

say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$CHAIN"; }
# pgrep (not ps|grep) so the pattern cannot match the matcher itself.
alive(){ pgrep -f "$1" >/dev/null 2>&1; }
# Waits while ANY of the given patterns matches. Watching only the per-run
# root.exe would be wrong: between two runs the xargs pool is momentarily empty
# and the stage would look finished, so the driver script is watched too.
wait_for(){ local what="$1"; shift; local n=0 any
  while :; do
    any=0; for p in "$@"; do alive "$p" && { any=1; break; }; done
    [ "$any" -eq 0 ] && break
    sleep 60; n=$((n+1))
    [ $((n % 30)) -eq 0 ] && say "  ... still waiting on $what (${n} min)"
  done
  say "$what finished"; }

say "=========== CHAIN 15C+p START ==========="

# ---- S1: FRIBDAQ reco --------------------------------------------------------
say "S1: waiting for FRIBDAQ reco"
wait_for "FRIBDAQ reco" "overnight_C15\.sh" "unpackReco_C15\.C"
say "S1 reco sizes:"
for r in $(cat "$HERE/runs_pp_fribdaq.txt"); do
  f="$RECO/${r}_reco.root"
  printf "    %-10s %s\n" "$r" "$([ -f "$f" ] && echo "$(( $(stat -c%s "$f") / 1048576 )) MB" || echo MISSING)" | tee -a "$CHAIN"
done

# ---- S2: first fit pass ------------------------------------------------------
say "S2: waiting for fit pass 1"
wait_for "fit pass 1" "fitall_C15\.sh" "fitUKF_C15\.C" "fitGenfit_C15\.C"

# ---- S3: fit the FRIBDAQ runs ------------------------------------------------
say "S3: fit pass 2 (FRIBDAQ runs)"
"$HERE/fitall_C15.sh" "" 6 >> "$LOG/fitall_pass2.log" 2>&1
say "S3 done: ukf=$(ls "$FIT"/*_ukf.root 2>/dev/null | wc -l) nomat=$(ls "$FIT"/*_genfit_nomat.root 2>/dev/null | wc -l) matFX=$(ls "$FIT"/*_genfit_mat.root 2>/dev/null | wc -l)"

# ---- S4: physics -------------------------------------------------------------
set +u
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"

# Runs that actually have a UKF fit, as a CSV for the analysis macros.
RUNSCSV=$(ls "$FIT"/*_ukf.root 2>/dev/null | sed 's#.*/##;s/_ukf.root//' | sort -u | paste -sd,)
say "S4: physics over $(echo $RUNSCSV | tr ',' ' ' | wc -w) runs"

# S4a: Ex spectrum per fitter at the placeholder Ebeam, so the three are comparable.
for fitter in ukf genfit_nomat genfit_mat; do
  say "S4a: Ex spectrum, fitter=$fitter (Ebeam=195 placeholder)"
  root -b -q -l "$HERE/pp/ex_C15.C(\"$RUNSCSV\",\"$FIT/\",195.0,5.0,\"_$fitter\",1.007825,15.0105993,\"\",\"$fitter\")" \
       > "$LOG/ex_${fitter}.log" 2>&1
done

# S4b: Ebeam scan on UKF. Ebeam=195 is a placeholder inherited from an old macro;
# the real value must put the elastic (p,p) peak at Ex = 0. Scan and read the peak
# position off each spectrum in the morning.
for E in 150 160 170 180 190 195 200 210; do
  say "S4b: Ebeam scan, UKF, Ebeam=$E"
  root -b -q -l "$HERE/pp/ex_C15.C(\"$RUNSCSV\",\"$FIT/\",${E}.0,5.0,\"_ukf_E${E}\",1.007825,15.0105993,\"\",\"ukf\")" \
       > "$LOG/ex_ukf_E${E}.log" 2>&1
done

# S4c: three-fitter comparison over the same reco.
say "S4c: compare fitters"
root -b -q -l "$HERE/pp/compare_fitters_C15.C(\"$RUNSCSV\",\"$FIT/\")" > "$LOG/compare_fitters.log" 2>&1

say "=========== CHAIN 15C+p DONE ==========="
say "plots in $HERE/pp/plots/ ; fits in $FIT"
