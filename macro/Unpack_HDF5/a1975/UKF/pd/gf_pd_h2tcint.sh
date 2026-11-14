#!/usr/bin/env bash
# =====================================================================================================
# 16C(p,d) genfit on the H2 FLOAT-TAU re-reco (reco_h2_tau/, HDBSCAN + PSAMultiFit float tau) with the
# deuteron gate on the TRACE-INTEGRAL plane (2026-11-13).
#
# ONE thing differs from the (p,d) fit in prod_h2_hdb.sh: the gate. Same fitter arguments verbatim
# (geometry ATTPC_H300torr_RT 3.308e-5, matEffects, CATIMA MSC+straggling, back-extrapolation,
# deuteron_H2_catima.txt, measSigma 4, pdg/mass), but
#   gate  = /mnt/f/a1975/reco_h2_tau/gate_pd_h2_integral.json (drawn by Yassid 2026-11-13 on
#           pid_plane_h2_integral.root, vertices 17-18 removed to undo a self-crossing)
#   gateTauNomTB = 4.5 -> AtGenfitter::SetPIDGateTauScaled: the gate sees AtSpyralPID's sqrtdEdx on
#           hit charge x tau_fit/4.5, Brho from the unscaled estimate -- the cache definition exactly.
# VERIFIED on run_0110 events 0-1999: the fitted set equals {cache in-gate} AND {fits with a
# pass-all gate} track for track (36 = 36, 0 extra, 0 missing).
#
# A run is marked .done only if the log ends with the fitter's "Done." line (a killed fit leaves a
# partial file and no segfault). The gate is copied next to the output for provenance.
#
#   NPAR=4 ./gf_pd_h2int.sh [runs...]
# =====================================================================================================
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/UKF
# config.sh reads unset variables: source it BEFORE set -u
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
export GENFIT=/home/yassid/fair_install/GenFit
export LD_LIBRARY_PATH=$GENFIT/lib:${LD_LIBRARY_PATH:-}
set -uo pipefail

export RECO=/mnt/f/a1975/reco_h2_tctau/
export OUT=/mnt/f/a1975/gf_pd_h2tcint/
export LOG=${OUT}logs/
export GATE=/mnt/f/a1975/reco_h2_tau/gate_pd_h2_integral.json
export GEO=ATTPC_H300torr_RT
export RHO=3.308e-5
export TAG=pdtcint
NPAR=${NPAR:-4}
mkdir -p "$OUT" "$LOG"
cp -p "$GATE" "${OUT}gate_used.json"

one() {
  run="run_$1"
  rc="${RECO}${run}_multifit_reco.root"
  fo="${OUT}${run}_multifit_genfitter_${TAG}.root"
  lg="${LOG}gf_${TAG}_${run}.log"
  [ -e "${fo}.done" ] && { echo "[have] $run"; return 0; }
  [ -e "${rc}.taudone" ] || { echo "[SKIP] $run: no finished float-tau reco"; return 1; }
  echo "[fit  ] $run $(date '+%H:%M:%S')"
  root -l -b -q "pipeline/fitGenfitter_a1975.C(\"${run}_multifit\",-1,\"$RECO\",\"_${TAG}\",\"$OUT\",\
-2.85,2,5,\"deuteron_H2_catima.txt\",\"$GATE\",4.0,kTRUE,1000010020,2.0135532,1,kFALSE,\"$GEO\",kTRUE,kTRUE,kFALSE,${RHO},kTRUE,\
kFALSE,kFALSE,4.5)" > "$lg" 2>&1
  if grep -qi 'segmentation violation' "$lg" || [ ! -s "$fo" ] || ! tail -3 "$lg" | grep -q 'Done\.'; then
    echo "[FAIL] $run (see $lg)"; return 1
  fi
  grep -q "PID gate on the TRACE-INTEGRAL plane" "$lg" || { echo "[FAIL] $run: integral-gate banner missing"; return 1; }
  grep -q "CATIMA material model: MSC ON, straggling ON" "$lg" || echo "[WARN] $run: CATIMA banner missing"
  grep -q "back-extrapolation to the beam axis: ON" "$lg" || echo "[WARN] $run: back-extrapolation banner missing"
  touch "${fo}.done"
  echo "[done ] $run $(date '+%H:%M:%S')  $(du -h "$fo" | cut -f1)"
}
export -f one

RUNS="${*:-$(seq -f %04g 106 189)}"
echo "=== (p,d) genfit, integral gate: $(echo $RUNS | wc -w) runs, $NPAR parallel, -> $OUT ==="
printf '%s\n' $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "=== gf_pd_h2int finished $(date '+%H:%M:%S') ==="
