#!/usr/bin/env bash
# Fit the a1975 16C(d,t)15C reconstruction at dv 1.10424 with BOTH fitters.
# The reco already exists (reco_d2_dv1104, 47/47); this is the fit stage only.
#
# TWO DELIBERATE DIFFERENCES FROM dv_dt.sh:
#
#   1. GATE.  pid/triton_d2_dv1104.json, drawn on the dv 1.10424 plane and on the IC-gated
#      sample. The old pid/triton_d2.json belongs to the dv 1.136 plane -- Brho goes as
#      1/|sin(polar)| and polar comes from a rho-vs-z regression, so a 2.9 % z rescale moves
#      every point under the polygon. A gate belongs to the plane it was drawn on.
#
#   2. GAS DENSITY 6.61e-5 FOR BOTH.  D2 at 300 torr, 293 K is 6.613e-5 g/cm3. The genfit
#      path already used that; fitUKF_a1975_deuterium.C DEFAULTS TO 9.0e-5, which is D2 at
#      about 408 torr. Left at the default the two fitters would be stopping tritons in
#      different gas, and the comparison would be meaningless.
#
# NOTE ON GATING, and it is not symmetric:
#   genfit applies the PID gate BEFORE fitting -- rejected tracks are never fitted.
#   the UKF macro applies NO gate; it fits everything and writes AtPIDEvent, so the gate is
#   applied at analysis time. The UKF output is therefore LARGER and a superset. Gate it in
#   the cache step, not here.
#
#   ./fit_dt_dv1104.sh [nparallel] [runs...]
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
# config.sh reads unset vars: source it BEFORE set -u or it kills the script silently
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export GENFIT=/home/yassid/fair_install/GenFit
export LD_LIBRARY_PATH=$GENFIT/lib:${LD_LIBRARY_PATH:-}
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
set -uo pipefail

NPAR="${1:-6}"; shift || true
RUNS="${*:-0016 0017 0018 0019 0020 0021 0022 0023 0026 0027 0031 0032 0034 0036 0037 0038 0039 0040 0041 0042 0043 0044 0046 0048 0057 0058 0076 0077 0078 0079 0080 0082 0083 0084 0085 0086 0087 0088 0089 0091 0092 0095 0096 0097 0098 0102 0103}"

export REC=/mnt/f/a1975/reco_d2_dv1104/
export GF=/mnt/f/a1975/gf_dt_dv1104/
export UK=/mnt/f/a1975/ukf_dt_dv1104/
export LOG=/mnt/f/a1975/logs_dv1104/
export PAR=ATTPC.a1975_deuterium_dv1104.par
export GATE=pid/triton_d2_dv1104.json
export RHO=6.61e-5
mkdir -p "$GF" "$UK" "$LOG"

one() {
  n="$1"; r="run_${n}"
  rc="${REC}${r}_multifit_reco.root"
  [ -s "$rc" ] || { echo "[noreco] $r"; return 0; }

  # ---- GENFIT (gate applied before the fit) --------------------------------------------
  fo="${GF}${r}_multifit_genfitter_t.root"
  if [ ! -s "${fo}.done" ]; then
    root -l -b -q "fitGenfitter_a1975_deuterium.C(\"${r}_multifit\",-1,\"$REC\",\"\",\"$GF\",-2.85,2,5,\"$GATE\",4.0,10.0,170.0,kFALSE,kTRUE,1000010030,3.01550072,1,\"t\",\"_reco\",\"ATTPC_D300torr_v2_geomanager.root\",kTRUE,${RHO},2,\"$PAR\")" \
      > "${LOG}gf_${r}.log" 2>&1
    if grep -q 'segmentation violation' "${LOG}gf_${r}.log" || [ ! -s "$fo" ]; then
      echo "[FAIL gf] $r"; rm -f "$fo"
    else
      touch "${fo}.done"; echo "[gf ok] $r"
    fi
  fi

  # ---- UKF (no gate; fits everything, writes AtPIDEvent) -------------------------------
  uo="${UK}${r}_genfitter_t_UKF.root"
  if [ ! -s "${uo}.done" ]; then
    root -l -b -q "fitUKF_a1975_deuterium.C(\"${r}\",-1,\"triton\",-1,2.85,${RHO},\"$REC\",\"$UK\",2.0,0.3,1,4,kTRUE,\"_multifit_reco\",\"_genfitter_t_UKF\")" \
      > "${LOG}ukf_${r}.log" 2>&1
    if grep -q 'segmentation violation' "${LOG}ukf_${r}.log" || [ ! -s "$uo" ]; then
      echo "[FAIL ukf] $r"; rm -f "$uo"
    else
      touch "${uo}.done"; echo "[ukf ok] $r"
    fi
  fi
  echo "[done] $r $(date '+%H:%M:%S')"
}
export -f one

printf '%s\n' $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "=== GENFIT $(ls ${GF}*_genfitter_t.root 2>/dev/null | wc -l)/47   UKF $(ls ${UK}*_genfitter_t_UKF.root 2>/dev/null | wc -l)/47 ==="
