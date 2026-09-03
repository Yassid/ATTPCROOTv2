#!/usr/bin/env bash
# ======================================================================================
# THE H2 (proton-target) RUNS ON THE STANDARD CONFIGURATION: AtPSAMultiFit + HDBSCAN.
# ======================================================================================
#
# WHY THIS EXISTS.  The a1975 H2 runs (0106-0189) feed three channels -- (p,p), (p,d) and (p,t) --
# and all three have been reconstructed since June with AtPSAMax + AtTrackFinderTC, because
# unpackReco_a1975_UKF.C hardcodes both and offers no option.  Meanwhile:
#
#   * the (d,t) DATA runs on AtPSAMultiFit + HDBSCAN (dv_dt.sh, mcs 0 / ms 3),
#   * the (p,d) and (p,t) SIMULATIONS run on AtPSAMultiFit + HDBSCAN (run_reco_C16pd.C,
#     run_reco_C16pt.C, mcs 20 / ms 8),
#
# so the two H2 transfer channels measure their yield with one track finder and their ACCEPTANCE
# with a different one.  That is not a bookkeeping detail: the looper study measured a 6.9x
# suppression of full-turn deuterons in the (p,d) data and none at all in (d,t), and the two
# candidate causes were exactly "different finder" and "different physics", never separated.
# This production removes the finder from the list.
#
# WHAT IS CHANGED, AND WHAT IS DELIBERATELY NOT.  One variable: PSA + pattern recognition.
#
#   PSA        AtPSAMax(thr 20)                 ->  AtPSAMultiFit(thr 20, peaking 0.720,
#                                                     maxPeaks 4, minSeparation 4)
#   cleaning   none                             ->  AtDirDeDxCleaner   (as sim, as (d,t) data)
#   PRA        AtTrackFinderTC(15.0, 7.5)       ->  AtTrackFinderHDBSCAN(mcs 20, ms 8, cse 10,
#                                                     mover, join 15, overlap 0.25, gap 40, ang 35)
#
#   space charge   STAYS ON.  It is a data-only distortion correction; the simulation has no
#                  space charge to undo, so "sim has none" is not a mismatch to fix.  The current
#                  H2 production has it on and removing it would move the vertex, which is a
#                  calibration change and not what is being tested here.
#   pad time offsets   STAY OFF.  pad_time_correction.csv was derived on the D2 runs and neither
#                  the current H2 production nor the H2 simulations use it.  Turning it on here
#                  would be a THIRD simultaneous change to the z of every hit.
#   parameter file     ATTPC.a1954.par, exactly as now.  Its DriftVelocity 1.30 / DriftLength
#                  1000 are what every existing H2 number is calibrated against, and they are
#                  identical to the sim par (ATTPC.a1975_C16_sim.par) -- the two par files differ
#                  only in Density and GasPressure, which the reco stage never reads.  So z is
#                  untouched and the comparison is clean.
#
#   mcs 20 / ms 8 and NOT the (d,t) data's mcs 0 / ms 3, because the point of the exercise is to
#   put the H2 DATA on the same finder configuration as the H2 ACCEPTANCE, and the H2 sims are
#   at 20/8.  (d,t) is internally consistent already and is not touched.
#
# NOTHING IS OVERWRITTEN.  New reco -> reco_h2_hdb/, new fits -> gf_pd_hdb/ and gf_pt_hdb/.  The
# TC production in reco/ and reco_pd_catima_bx/ and reco_pt_catima/ stays exactly where it is, so
# the two can be differenced run by run.  That difference IS the deliverable.
#
# The genfit call is copied verbatim from prod_pd_catima.sh / prod_pt_catima.sh -- same geometry
# (ATTPC_H300torr_RT, 3.308e-5), same CATIMA MSC+straggling, same back-extrapolation, same gate,
# same dE/dx table.  Only the input reco directory and the output suffix change.  If the fitter
# arguments here ever drift from those two files the comparison stops meaning anything.
#
#   ./prod_h2_hdb.sh [nparallel] [runs...]
#
# Runs are ordered so the six runs the hit-geometry studies use (0106-0110, 0113) complete first
# and the consistency check can start while the rest are still going.
# ======================================================================================
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/UKF
# config.sh reads unset vars: source it BEFORE set -u or it kills the script silently
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export GENFIT=/home/yassid/fair_install/GenFit
export LD_LIBRARY_PATH=$GENFIT/lib:${LD_LIBRARY_PATH:-}
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
set -uo pipefail

NPAR="${1:-5}"; shift || true

export H5=/mnt/f/a1975/h5/
export RECO=/mnt/f/a1975/reco_h2_hdb/
export PDOUT=/mnt/f/a1975/gf_pd_hdb/
export PTOUT=/mnt/f/a1975/gf_pt_hdb/
export LOG=/mnt/f/a1975/logs_h2_hdb/
export PAR=ATTPC.a1954.par
export GEO=ATTPC_H300torr_RT
export RHO=3.308e-5
mkdir -p "$RECO" "$PDOUT" "$PTOUT" "$LOG"

one() {
  r="$1"; run="run_${r}"
  rc="${RECO}${run}_multifit_reco.root"

  # ---------------- stage 1: reco (unpack + MultiFit + SC + clean + HDBSCAN) -------------
  # A marker means COMPLETED, never "a file exists": the macro's own Done line, a size floor and
  # no segfault.  A failed stage deletes its own partial output so a rerun redoes it.
  # -e, NOT -s: the marker is created by `touch`, so it is ZERO BYTES and -s is false for it.
  # With -s this guard never fires and a rerun silently redoes every completed run from scratch.
  # (dv_dt.sh carries the same latent bug.)
  if [ ! -e "${rc}.done" ]; then
    # persistRaw = FALSE.  AtRawEvent is 10x the rest of the file: with it on, one run projected to
    # ~16 GB and the 84-run set to 1.3 TB, against ~1.2 GB per run without.  Nothing downstream
    # reads it -- the fitter takes AtPatternEvent, the ion chamber comes from the separate FRIB
    # file -- and AtEventH is persisted unconditionally by the macro, so the full PSA hit cloud
    # (what the topology study needs) is kept either way.  dv_dt.sh passes false for the same reason.
    root -l -b -q "../D2_UKF/unpackReco_multifit.C(\"$run\",0,false,\"$RECO\",\"$H5\",true,false,\"multifit\",0,20,\"hdbscan\",20,8,0,0.1,\"$PAR\",kTRUE)" \
        > "${LOG}reco_${run}.log" 2>&1
    if grep -qi 'segmentation violation' "${LOG}reco_${run}.log" || [ ! -s "$rc" ] \
       || [ "$(stat -c%s "$rc" 2>/dev/null || echo 0)" -lt 20000000 ]; then
      echo "[FAIL reco] $run  (see ${LOG}reco_${run}.log)"; rm -f "$rc"; return 1
    fi
    # the banner is the only proof the finder actually was HDBSCAN and not the macro's tc default
    grep -q "AtTrackFinderHDBSCAN mover (mcs 20, ms 8" "${LOG}reco_${run}.log" \
      || { echo "[FAIL reco] $run: HDBSCAN banner missing -- wrong finder"; rm -f "$rc"; return 1; }
    grep -q "PSA   : AtPSAMultiFit" "${LOG}reco_${run}.log" \
      || { echo "[FAIL reco] $run: MultiFit banner missing"; rm -f "$rc"; return 1; }
    touch "${rc}.done"; echo "[reco ok] $run  $(date '+%H:%M:%S')"
  fi

  # ---------------- stage 2: genfit, (p,d) then (p,t) -----------------------------------
  # Arguments verbatim from prod_pd_catima.sh / prod_pt_catima.sh.  fitGenfitter_a1975.C builds
  # its input name as <recoDir><run>_reco.root, but this reco is named <run>_multifit_reco.root,
  # so the run string passed to the fitter carries the _multifit tag.
  # tail carries the per-channel trailing arguments verbatim from that channel's own driver:
  # (p,t) was migrated to native CATIMA dE/dx on 2026-08-19 and (p,d) never was, and this is not
  # the place to close that gap -- one variable at a time.
  fitone() {
    out="$1"; sfx="$2"; tab="$3"; gate="$4"; pdg="$5"; mass="$6"; tag="$7"; tail="$8"
    fo="${out}${run}_multifit_genfitter_${tag}.root"
    [ -e "${fo}.done" ] && return 0   # -e: the marker is a zero-byte touch file
    root -l -b -q "pipeline/fitGenfitter_a1975.C(\"${run}_multifit\",-1,\"$RECO\",\"${sfx}\",\"$out\",\
-2.85,2,5,\"$tab\",\"$gate\",4.0,kTRUE,${pdg},${mass},1,kFALSE,\"$GEO\",kTRUE,kTRUE,kFALSE,${RHO},kTRUE${tail})" \
      > "${LOG}gf_${tag}_${run}.log" 2>&1
    if grep -qi 'segmentation violation' "${LOG}gf_${tag}_${run}.log" || [ ! -s "$fo" ]; then
      echo "[FAIL fit $tag] $run  (see ${LOG}gf_${tag}_${run}.log)"; rm -f "$fo"; return 1
    fi
    grep -q "CATIMA material model: MSC ON, straggling ON" "${LOG}gf_${tag}_${run}.log" \
      || echo "[WARN] $run $tag: CATIMA banner missing"
    grep -q "back-extrapolation to the beam axis: ON" "${LOG}gf_${tag}_${run}.log" \
      || echo "[WARN] $run $tag: back-extrapolation banner missing"
    touch "${fo}.done"; echo "[fit $tag ok] $run  $(date '+%H:%M:%S')"
  }

  # NOTE the (p,t) ejectile mass: 3.01604928 is the ATOMIC mass of tritium, while (d,t) fits with
  # 3.01550072, the NUCLEAR (triton) mass -- a 0.00055 u = 0.51 MeV difference between the two
  # channels' tritons.  That inconsistency is inherited verbatim from prod_pt_catima.sh ON PURPOSE:
  # this production exists to isolate the track finder, and changing the mass at the same time
  # would make the difference unattributable.  Flagged, not fixed here.
  fitone "$PDOUT" "_pdhdb" "deuteron_H2_catima.txt" "pid/deuteron_band.json" 1000010020 2.0135532 "pdhdb" ""
  fitone "$PTOUT" "_pthdb" "triton_H2_catima.txt" \
         "/home/yassid/a1975_C16dt_analysis/triton_pt.json" 1000010030 3.01604928 "pthdb" ",kTRUE,kFALSE"
}
export -f one

# The six runs the hit-geometry and winding studies use, first; then the rest in order.
FIRST="0106 0107 0108 0109 0110 0113"
REST=$(seq -f "%04g" 106 189 | grep -vwE "$(echo $FIRST | tr ' ' '|')")
RUNS="${*:-$FIRST $REST}"

echo "=== H2 standard-configuration production: $(echo $RUNS | wc -w) runs, $NPAR parallel ==="
echo "=== PSA AtPSAMultiFit(20) + AtDirDeDxCleaner + HDBSCAN(mcs 20, ms 8), SC ON, par $PAR ==="
echo "=== reco -> $RECO   (p,d) -> $PDOUT   (p,t) -> $PTOUT ==="
printf '%s\n' $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "=== done: $(ls "$RECO"/*_multifit_reco.root 2>/dev/null | wc -l) reco, \
$(ls "$PDOUT"/*.root 2>/dev/null | wc -l) pd fits, $(ls "$PTOUT"/*.root 2>/dev/null | wc -l) pt fits ==="
