#!/usr/bin/env bash
# Build the (d,t) kinematics cache from the UKF production (ukf_dt_full).
#
# Same gate, same Ebeam, same "no cuts" policy as cache_dt_dv1104.sh and cache_dt_maton.sh, so
# the three caches differ ONLY by the fitter and can be dropped into one explorer page.
#
# TWO THINGS THAT ARE DIFFERENT ABOUT THE UKF AND BOTH BITE SILENTLY:
#
#   1. FILE NAMING.  The UKF writes run_XXXX_genfitter_t_UKF.root -- no "_multifit" stage tag and
#      a "_UKF" suffix -- whereas genfit writes run_XXXX_multifit_genfitter_t.root. ex_dt_a1975.C
#      builds the path as inDir + run + "_genfitter_t" + suffix, so the run string here is plain
#      run_XXXX and the suffix is "_UKF".
#
#   2. THE KINEMATICS SLOTS ARE SWAPPED between the two fitters:
#           AtGenfitter   SetKinematics(KEfit)=RAW at the first measurement point
#                         SetKinematicsXtr(KE) =back-extrapolated to the beam axis  <- corrected
#           AtFitterUKF   SetKinematics(KE)    =at the vertex                       <- corrected
#                         SetKinematicsXtr(...)=raw, at the first cluster
#      ex_dt_a1975.C writes Xtr into `ke`/`theta`/`ex` and the other slot into
#      `kefit`/`thetafit`/`exfit`. So for GENFIT the corrected columns are ke/theta, and for the
#      UKF they are kefit/thetafit. This cache is written either way -- BOTH column pairs are
#      present -- but anything reading it must take kefit/thetafit for the UKF. The explorer
#      script does exactly that via mkexp_pp's keCol/thCol arguments.
#
# The UKF macro applies NO PID gate of its own; the gate is applied here, which is what makes
# this sample comparable to the genfit ones (genfit gates before fitting).
#
#   ./cache_dt_ukf.sh [nparallel] [workdir]
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
NPAR="${1:-6}"
WORK="${2:-/mnt/f/a1975/caches/.ukf}"
OUT=/mnt/f/a1975/caches/dt_kin_ukf.root
FIT=/mnt/f/a1975/ukf_dt_full/
EBEAM=184.17
mkdir -p "$WORK" /mnt/f/a1975/caches
# config.sh reads unset vars; under `set -u` it kills the script SILENTLY, zero-byte log.
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
set -uo pipefail
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}

NUMS="0016 0017 0018 0019 0020 0021 0022 0023 0026 0027 0031 0032 0034 0036 0037 0038 0039 0040 0041 0042 0043 0044 0046 0048 0057 0058 0076 0077 0078 0079 0080 0082 0083 0084 0085 0086 0087 0088 0089 0091 0092 0095 0096 0097 0098 0102 0103"

one() {
  n="$1"; work="$2"; fit="$3"; eb="$4"
  out="${work}/run_${n}.root"; log="${work}/log_${n}.txt"
  [ -s "$out" ] && { echo "[skip] run_$n"; return 0; }
  [ -s "${fit}run_${n}_genfitter_t_UKF.root" ] || { echo "[nofit] run_$n"; return 0; }
  root -b -l -q "ex_dt_a1975.C(\"run_${n}\",\"${fit}\",\"_UKF\",\"pid/triton_d2_dv1104.json\",${eb},1e9,0.0,1e9,0.0,180.0,\"${out}\",\"${work}/plot_${n}.png\")" \
     > "$log" 2>&1
  if grep -qi "no IC file" "$log"; then echo "[NO-IC] run_$n  <-- beam gate missing"; fi
  [ -s "$out" ] && echo "[done] run_$n $(date '+%H:%M:%S')" || echo "[FAIL] run_$n"
}
export -f one

printf '%s\n' $NUMS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {} "$WORK" "$FIT" "$EBEAM"

n=$(ls "$WORK"/run_*.root 2>/dev/null | wc -l)
echo "=== $n per-run caches; merging ==="
hadd -f "$OUT" "$WORK"/run_*.root > "$WORK/hadd.log" 2>&1 && echo "=== cache: $OUT ===" || echo "=== HADD FAILED, see $WORK/hadd.log ==="
