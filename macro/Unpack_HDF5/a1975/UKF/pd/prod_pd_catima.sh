#!/usr/bin/env bash
# 16C(p,d)15C production WITH material effects, using the CATIMA backend.
#
# This configuration only became usable on 2026-08-17. genfit's TGeoMaterialInterface navigates
# the one shared gGeoManager and never restores its state, so a single track whose extrapolation
# wandered left the navigator parked and poisoned every LATER fit in the job -- on the (d,t) set
# that cost 58.6% of all tracks. AtGenfitter now re-seeds the navigator before each fit.
#
# BACK-EXTRAPOLATION, added 2026-08-17 and the reason for the _bx output.
#   The first CATIMA production (reco_pd_catima/, _pdcat) ran with NO back-extrapolation:
#   fitGenfitter_a1975.C had no such argument and AtGenfitter defaults fBackExtrapToAxis to
#   kFALSE, so GetKinematicsXtr() came back equal to the raw fit for 456/456 tracks of run_0106.
#   The vertex gap -- the gas the deuteron crossed before its first cluster -- was therefore not
#   corrected at all. The last argument is now backExtrap = kTRUE, which lets genfit transport
#   the state to the beam axis through CATIMA material.
#
#   The gap is corrected ONCE, by genfit. There is no manual dE/dx term here and there must not
#   be: AtGenfitter's SetManualELoss block runs unconditionally, so pairing it with matEffects
#   counts the gap twice. On (d,t) that inflated the low-branch correction from +4.65% to +9.5%.
#   fitGenfitter_a1975.C exposes no manual knob, which is what keeps this safe.
#
#   Only GetKinematicsXtr() changes. The raw fit is bit-identical to the _pdcat production, so
#   the two can be differenced to measure the gap correction directly.
#
# NOTHING IS OVERWRITTEN. Output goes to reco_pd_catima_bx/ with suffix _pdcatbx, alongside both
# the first CATIMA production (reco_pd_catima/, _pdcat) and the matFX-off one in reco_pd/ (_pd).
#
# THE GEOMETRY IS THE THING TO GET RIGHT HERE.
#   fitGenfitter_a1975.C defaults to geoName = ATTPC_H1bar, which is 8.27e-5 g/cm3. The a1975
#   gas is H2 at 300 torr = 3.308e-5 (ATTPC_H300torr_RT). That default was HARMLESS while
#   matEffects was off, because with no material effects the density cannot reach the fit --
#   which is exactly why it survived unnoticed. Turning material effects on makes it live, and
#   8.27e-5 is 2.5x too much material. Verified from the geometry files themselves:
#     ATTPC_H1bar        8.270e-5 g/cm3   <- wrong for a1975
#     ATTPC_H300torr     3.553e-5 g/cm3   <- also not it
#     ATTPC_H300torr_RT  3.308e-5 g/cm3   <- correct
#
# OTHER SETTINGS AND WHY
#   catimaMSC / catimaStraggling  BOTH on. On (d,t) the collapse was a conjunction: removing
#       either one, or the dE/dx table, put it back to ~50%.
#   elossDensity 3.308e-5   loads the dE/dx table in HYBRID mode. Below beta*gamma = 0.05
#       (KE = 2.3 MeV for a deuteron) genfit applies NO stopping power without a curve loaded.
#   matFallback kFALSE      a throwing track is dropped, not silently refitted with
#       setNoEffects(true) and kept, which would blend two physics models into one spectrum.
#
#   ./prod_pd_catima.sh [nparallel]
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/UKF
# config.sh reads unset vars: source it BEFORE set -u or it kills the script silently
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export GENFIT=/home/yassid/fair_install/GenFit
export LD_LIBRARY_PATH=$GENFIT/lib:${LD_LIBRARY_PATH:-}
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
set -uo pipefail

NPAR="${1:-4}"                      # genfit fits are CPU-bound at ~0.8 GB RSS each
export RECO=/mnt/f/a1975/reco/
export OUT=/mnt/f/a1975/reco_pd_catima_bx/
export LOG=/mnt/f/a1975/logs_pd_catima_bx/
export GEO=ATTPC_H300torr_RT
export RHO=3.308e-5
export TAB=deuteron_H2_catima.txt
export GATE=pid/deuteron_band.json
mkdir -p "$OUT" "$LOG"

one() {
  r="$1"; run="run_${r}"
  [ -s "${RECO}${run}_reco.root" ] || { echo "[noreco] $run"; return 0; }
  fo="${OUT}${run}_genfitter_pdcatbx.root"
  [ -s "${fo}.done" ] && { echo "[have] $run"; return 0; }
  root -l -b -q "pipeline/fitGenfitter_a1975.C(\"${run}\",-1,\"$RECO\",\"_pdcatbx\",\"$OUT\",\
-2.85,2,5,\"$TAB\",\"$GATE\",4.0,kTRUE,1000010020,2.0135532,1,kFALSE,\"$GEO\",kTRUE,kTRUE,kFALSE,${RHO},kTRUE)" \
    > "${LOG}gf_${run}.log" 2>&1
  if grep -qi 'segmentation violation' "${LOG}gf_${run}.log" || [ ! -s "$fo" ]; then
    echo "[FAIL] $run  (see ${LOG}gf_${run}.log)"; rm -f "$fo"
  else
    # catch a run whose backend or geometry silently did not engage, here and not at analysis time
    grep -q "CATIMA material model: MSC ON, straggling ON" "${LOG}gf_${run}.log" \
      || echo "[WARN] $run: CATIMA banner missing"
    grep -qi "WARNING: matEffects with geoName" "${LOG}gf_${run}.log" \
      && echo "[WARN] $run: wrong geometry for material effects"
    # the whole point of this production: a missing banner means it silently ran as _pdcat again
    grep -q "back-extrapolation to the beam axis: ON" "${LOG}gf_${run}.log" \
      || echo "[WARN] $run: back-extrapolation banner missing"
    touch "${fo}.done"; echo "[ok] $run  $(date '+%H:%M:%S')"
  fi
}
export -f one

# seq -w pads to the WIDTH OF THE LARGEST ARGUMENT, i.e. 3 digits -> "106", while the reco
# files are run_0106_reco.root. genfit_production.sh has the same latent bug.
RUNS=$(seq -f "%04g" 106 189)
echo "=== (p,d) CATIMA production: $(echo $RUNS | wc -w) runs, $NPAR parallel -> $OUT ==="
echo "=== geometry $GEO (rho $RHO g/cm3), matFX ON, CATIMA MSC+straggling, fallback OFF ==="
printf '%s\n' $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "=== done: $(ls "$OUT"/*_genfitter_pdcatbx.root 2>/dev/null | wc -l) fit files in $OUT ==="
