#!/usr/bin/env bash
# =====================================================================================================
# 16C(d,t) genfit on the D2 HDBSCAN + FLOAT-TAU re-reco (reco_d2_hdbtau/) with the triton gate on the
# TRACE-INTEGRAL plane (2026-11-18).
#
# Fitter arguments = the (d,t) r280 production (sim_rebuild/r280/run_prod_dt.sh: CATIMA MSC+straggling,
# hybrid dE/dx table triton_D2_300torr.txt, matEffects, back-extrapolation, measSigma 4, -2.85 T,
# r280 geometry, par dv1104), except
#   reco  = /mnt/f/a1975/reco_d2_hdbtau/
#   gate  = gate_dt_d2_integral.json, drawn by Yassid 2026-11-18 on pid_plane_d2_integral.root
#           (frozen copy: ${OUT}gate_used.json -- the fit reads the COPY, so a re-save in the drawer cannot change it)
#   gateTauNomTB = 4.5 -> AtGenfitter::SetPIDGateTauScaled, the pid_width_dp.C / cache definition.
# VERIFIED run_0018 entries 0-1999: 60 gated fits (log) = 60 unique fitted KE in the tree, all 60 found
# with identical KE in a pass-all-gate fit of the same events; 82 plane tracks in the gate, the other 22
# fail after the gate (fit/theta/vertex). Tree entries are NOT reco entries (duplicate-entry trap).
#
#   NPAR=6 ./gf_dt_hdbtau_int.sh [runs...]      resume = rerun
# =====================================================================================================
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
set -uo pipefail

export RECO=/mnt/f/a1975/reco_d2_hdbtau/
export OUT=/mnt/f/a1975/gf_dt_hdbtau_int/
export LOG=${OUT}logs/
export GATE=${OUT}gate_used.json
NPAR=${NPAR:-6}
mkdir -p "$OUT" "$LOG"
exec 9>"${OUT}.lock"
flock -n 9 || { echo "gf_dt_hdbtau_int.sh: another instance holds the lock -- abort"; exit 1; }
[ -s "$GATE" ] || cp -p /mnt/f/a1975/reco_d2_hdbtau/gate_dt_d2_integral.json "$GATE"

one() {
  r="run_$1"
  fo="${OUT}${r}_multifit_genfitter_t.root"
  lg="${LOG}gf_${r}.log"
  if [ -e "${fo}.done" ] && [ -s "$fo" ]; then echo "[have] $r"; return 0; fi
  [ -e "${RECO}${r}_multifit_reco.root.taudone" ] || [ -s "${RECO}${r}_multifit_reco.root" ] || { echo "[SKIP] $r: no reco"; return 1; }
  echo "[fit ] $r $(date '+%H:%M:%S')"
  root -l -b -q "fitGenfitter_a1975_deuterium.C(\"${r}_multifit\",-1,\"$RECO\",\"\",\"$OUT\",-2.85,2,5,\"$GATE\",4.0,10.0,170.0,kTRUE,kTRUE,1000010030,3.01550072,1,\"t\",\"_reco\",\"ATTPC_D300torr_r280_geomanager.root\",kTRUE,0,2,\"ATTPC.a1975_deuterium_dv1104.par\",kFALSE,kFALSE,kFALSE,\"triton_D2_300torr.txt\",kTRUE,kTRUE,kTRUE,kFALSE,0,kFALSE,0,kFALSE,5.0,25.0,0,95.0,8,kFALSE,kFALSE,\"\",kFALSE,4.5)" \
    > "$lg" 2>&1 < /dev/null
  vr=$(root -b -l -q "pid/root_ok.C(\"$fo\")" 2>/dev/null | grep -E '^(VALID|INVALID)' | head -1)
  if [[ "$vr" == VALID* ]] && grep -q "Done\." "$lg" && grep -q "dE/dx from CATIMA" "$lg" && grep -q "TRACE-INTEGRAL" "$lg" \
     && ! grep -qi 'segmentation violation' "$lg"; then
    touch "${fo}.done"; echo "[ok  ] $r $(date '+%H:%M:%S') $vr"
  else
    echo "[FAIL] $r ${vr:-unreadable} (see $lg)"; rm -f "$fo"
  fi
}
export -f one

if [ $# -gt 0 ]; then printf '%s\n' "$@"; else
  ls ${RECO}run_*_multifit_reco.root | sed -E 's#.*/run_([0-9]{4})_multifit_reco.root#\1#'; fi | xargs -P "$NPAR" -I{} bash -c 'one {}'
echo "=== gf_dt_hdbtau_int finished $(date '+%F %T'): $(ls ${OUT}*.done 2>/dev/null | wc -l) done ==="
