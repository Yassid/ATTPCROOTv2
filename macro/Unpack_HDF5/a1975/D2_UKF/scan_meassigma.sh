#!/usr/bin/env bash
# measSigma scan for the (d,t) low branch.
#
# HYPOTHESIS: the ~2 MeV floor in reconstructed KE comes from the fit, not from the selection.
# AtGenfitter builds the per-cluster covariance as varT = sigma^2 + D_T^2 * L_drift, so
# measSigma is the FLOOR of the measurement uncertainty. A g.s. triton at 15 deg has a
# cyclotron radius of 81 mm; 4 mm on such an arc is ~5%, and with few clusters the curvature
# is then so weakly constrained that the fit leans on its seed -- which biases the momentum
# high and flattens the KE-vs-theta locus. In (p,d) the same 4 mm sat on a 312 mm radius and
# was irrelevant, which is why that channel never showed this.
#
# Everything except measSigma is identical to the production fit.
#   ./scan_meassigma.sh [nparallel] [runs...]
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export GENFIT=/home/yassid/fair_install/GenFit
export LD_LIBRARY_PATH=$GENFIT/lib:${LD_LIBRARY_PATH:-}
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
set -uo pipefail

NPAR="${1:-3}"; shift || true
RUNS="${*:-0016 0019 0022}"
export REC=/mnt/f/a1975/reco_d2_dv1104/
export PAR=ATTPC.a1975_deuterium_dv1104.par
export GATE=pid/triton_d2_dv1104.json
export LOG=/mnt/f/a1975/logs_ms/
mkdir -p "$LOG"

one() {
  spec="$1"; ms="${spec%%:*}"; n="${spec##*:}"; r="run_${n}"
  tag=$(echo "$ms" | tr -d '.')
  out=/mnt/f/a1975/gf_dt_ms${tag}/; mkdir -p "$out"
  fo="${out}${r}_multifit_genfitter_t.root"
  [ -s "${fo}.done" ] && { echo "[skip] ms=$ms $r"; return 0; }
  root -l -b -q "fitGenfitter_a1975_deuterium.C(\"${r}_multifit\",-1,\"$REC\",\"\",\"$out\",-2.85,2,5,\"$GATE\",${ms},10.0,170.0,kFALSE,kTRUE,1000010030,3.01550072,1,\"t\",\"_reco\",\"ATTPC_D300torr_v2_geomanager.root\",kTRUE,6.61e-5,2,\"$PAR\")" \
    > "${LOG}ms${tag}_${r}.log" 2>&1
  if grep -q 'segmentation violation' "${LOG}ms${tag}_${r}.log" || [ ! -s "$fo" ]; then
    echo "[FAIL] ms=$ms $r"
  else touch "${fo}.done"; echo "[ok] ms=$ms $r $(date '+%H:%M:%S')"; fi
}
export -f one

JOBS=""
for ms in 4.0 2.0 1.0; do for n in $RUNS; do JOBS="$JOBS ${ms}:${n}"; done; done
printf '%s\n' $JOBS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "=== scan done ==="
for ms in 40 20 10; do echo "  ms$ms: $(ls /mnt/f/a1975/gf_dt_ms${ms}/*_genfitter_t.root 2>/dev/null|wc -l) fits"; done
