#!/bin/bash
# ============================================================================
# a2091 15C(p,p') overnight #2 -- 2026-07-29. Runs on the IC + PID GATED data.
#
#   S1  wait for the gated pass (gatefit_C15.sh) to finish
#   S2  re-calibrate Ebeam on the GATED sample (pp/ebeam_C15.C, ridge + tilt)
#   S3  GENFIT noMat + matFX on the GATED recos -- cheap now, the gated files are small
#   S4  Ex spectra for all three fitters at the re-calibrated Ebeam + fitter comparison
#   S5  write a plain-text report
#
# Everything is idempotent (skips existing outputs) and writes to the Seagate-symlinked
# fit directory. Log lives OUTSIDE any directory that might get migrated.
# ============================================================================
set -u
REPO=/home/yassid/fair_install/ATTPCROOTv2
HERE=$REPO/macro/Unpack_HDF5/a2091/UKF
GIN=/home/yassid/a2091_C15_fit/in
GFIT=/home/yassid/a2091_C15_fit/gated
LOG=/home/yassid/a2091_C15_fit/logs/gated
MASTER=/home/yassid/a2091_C15_fit/overnight2.log
GEOM=ATTPC_H300torr_RT
DENS=3.308e-5
NPAR=6
mkdir -p "$GFIT" "$LOG"

say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$MASTER"; }
alive(){ pgrep -f "$1" >/dev/null 2>&1; }     # pgrep, not ps|grep: grep matches its own cmdline
trap 'pkill -P $$ 2>/dev/null' EXIT INT TERM  # killing me must kill the xargs pool too

# set +u around these: sourcing thisroot.sh under `set -u` hits an unbound variable and EXITS
srcbuild(){ set +u
            source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
            source "$REPO/$1/config.sh" >/dev/null 2>&1
            export ROOT_INCLUDE_PATH="$REPO/$1/include:/home/yassid/fair_install/FairRootInstall/include"
            set -u; }

say "=========== OVERNIGHT#2 (gated) START ==========="

# ---- S1: wait for the gated pass -------------------------------------------------
say "S1: waiting for the gated pass"
n=0
while alive "gatefit_C15\.sh" || alive "gate_events_C15\.C" || alive "fitUKF_C15\.C"; do
  sleep 60; n=$((n+1)); [ $((n % 30)) -eq 0 ] && say "  ... still gating (${n} min)"
done
say "S1 done: $(ls "$GIN"/*_reco.root 2>/dev/null | wc -l) gated recos, $(ls "$GFIT"/*_ukf.root 2>/dev/null | wc -l) UKF fits"

RUNS=$(ls "$GIN"/*_reco.root 2>/dev/null | sed 's#.*/##;s/_reco.root//' | sort -u | tr '\n' ' ')
[ -z "$RUNS" ] && { say "FATAL: no gated recos, nothing to do"; exit 1; }

# ---- S2: Ebeam on the gated sample ----------------------------------------------
srcbuild build
CSV=$(ls "$GFIT"/*_ukf.root 2>/dev/null | sed 's#.*/##;s/_ukf.root//' | sort -u | paste -sd,)
say "S2: Ebeam calibration on the gated sample"
# ALWAYS regenerate the gated proton_kin cache over the full run set. Do not reuse an existing
# one: a 2-run smoke-test cache was sitting there, and ebeam_C15.C would have silently
# calibrated on 353 tracks instead of the whole gated sample.
root -b -q -l "$HERE/pp/ex_C15.C(\"$CSV\",\"$GFIT/\",157.0,5.0,\"_gated\",1.007825,15.0105993,\"gated\",\"ukf\")" \
     > "$LOG/ex_gated_seed.log" 2>&1
root -b -q -l "$HERE/pp/ebeam_C15.C(\"plots/proton_kin_gated.root\")" > "$LOG/ebeam_gated.log" 2>&1
EB=$(grep -oP 'CONSENSUS Ebeam = \K[0-9]+' "$LOG/ebeam_gated.log" | tail -1)
if [ -z "${EB:-}" ]; then EB=157; say "S2: no consensus parsed, falling back to Ebeam=157"; fi
say "S2 done: Ebeam = $EB MeV"
grep -E "RIDGE:|TILT:|CONSENSUS|WARNING|tracks:" "$LOG/ebeam_gated.log" | tee -a "$MASTER"

# ---- S3: GENFIT on the gated recos ----------------------------------------------
srcbuild build_genfit
gf(){ local r="$1" tag="$2" mat="$3"
  [ -f "$GFIT/${r}_genfit${tag}.root" ] && return
  root -b -q -l "$HERE/pipeline/fitGenfit_C15.C(\"$r\",-1,\"$GIN/\",\"$tag\",\"$GFIT/\",-2.85,2,5,\"\",4.0,10.0,170.0,$mat,kFALSE,\"proton\",\"$GEOM\")" \
       > "$LOG/${r}_genfit${tag}.log" 2>&1
}
export -f gf; export HERE GIN GFIT LOG GEOM
say "S3: GENFIT matEffects=OFF on gated recos"
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'gf "$1" _nomat kFALSE' _ {} >> "$MASTER" 2>&1
say "S3a done: $(ls "$GFIT"/*_genfit_nomat.root 2>/dev/null | wc -l) files"
say "S3: GENFIT matEffects=ON on gated recos"
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'gf "$1" _mat kTRUE' _ {} >> "$MASTER" 2>&1
say "S3b done: $(ls "$GFIT"/*_genfit_mat.root 2>/dev/null | wc -l) files"

# ---- S4: Ex per fitter at the calibrated Ebeam ----------------------------------
srcbuild build
for fitter in ukf genfit_nomat genfit_mat; do
  C=$(ls "$GFIT"/*_${fitter}.root 2>/dev/null | sed "s#.*/##;s/_${fitter}.root//" | sort -u | paste -sd,)
  [ -z "$C" ] && { say "S4: no $fitter fits, skipping"; continue; }
  say "S4: Ex spectrum, $fitter, Ebeam=$EB"
  root -b -q -l "$HERE/pp/ex_C15.C(\"$C\",\"$GFIT/\",${EB}.0,5.0,\"_g_$fitter\",1.007825,15.0105993,\"15C(p,p') gated\",\"$fitter\")" \
       > "$LOG/ex_g_${fitter}.log" 2>&1
  grep -E "protons -> Ex" "$LOG/ex_g_${fitter}.log" | tail -1 | tee -a "$MASTER"
done
say "S4: fitter comparison"
root -b -q -l "$HERE/pp/compare_fitters_C15.C(\"$(echo $RUNS | tr ' ' ',')\",\"$GFIT/\")" > "$LOG/compare_gated.log" 2>&1

# ---- S5: report -----------------------------------------------------------------
REP=/home/yassid/a2091_C15_fit/OVERNIGHT2_REPORT.txt
{
  echo "a2091 15C(p,p') -- gated production, $(date)"
  echo
  echo "GATES   IC (15C beam) : pid/ic_15C.json  979.5 - 1278.8 ADC, single pulse"
  echo "        proton PID    : pid/proton_15C.json (19 vertices, NO arclen cut)"
  echo "        polar         : > 90 deg (PID convention, = 180 - theta_fitted)"
  echo "EBEAM   $EB MeV  ($(python3 -c "print(f'{$EB/15:.2f}')") MeV/u)  -- re-calibrated on the gated sample"
  echo
  echo "COUNTS"
  printf "  gated recos      : %s\n" "$(ls "$GIN"/*_reco.root 2>/dev/null | wc -l)"
  for t in ukf genfit_nomat genfit_mat; do printf "  %-16s : %s\n" "$t" "$(ls "$GFIT"/*_${t}.root 2>/dev/null | wc -l)"; done
  echo
  echo "EBEAM CALIBRATION"; grep -E "RIDGE:|TILT:|CONSENSUS|WARNING" "$LOG/ebeam_gated.log" 2>/dev/null | sed 's/^/  /'
  echo
  echo "EX YIELDS"
  for t in ukf genfit_nomat genfit_mat; do
    printf "  %-16s : %s\n" "$t" "$(grep -oE 'protons -> Ex: [0-9]+' "$LOG/ex_g_${t}.log" 2>/dev/null | tail -1)"
  done
  echo
  echo "PLOTS  $HERE/pp/plots/ex_C15_g_*.png , compare_fitters_C15.png"
  echo "       $HERE/pid/plots/ (PID plane, IC spectrum, gates)"
  echo
  echo "OPEN   - 15C has ONE bound excited state, 5/2+ at 0.740 MeV (Sn = 1.218):"
  echo "         check whether it now separates from the elastic tail in the gated Ex spectrum."
  echo "       - the ungated data had a KE ~ 0.79 MeV reconstruction threshold and a polar"
  echo "         88-94 deg blind spot; verify the gated theta window avoids both."
  echo "       - 39 UNGATED ukf fits + 23 genfit_nomat remain on the Seagate for comparison."
} > "$REP"
say "S5: report -> $REP"
say "=========== OVERNIGHT#2 DONE ==========="
