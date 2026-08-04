#!/usr/bin/env bash
# Build the (d,t) kinematics cache for the dv = 1.35 production.
#
# Ebeam is 184.0 MeV = 11.5 MeV/u, the correct beam energy. The older caches were built at 180
# and the free-Ebeam fits kept running to ~204, which is above the nominal incident energy and
# therefore an artefact (a momentum/B-field scale error being absorbed), not a beam energy.
#
# The cache carries BOTH energies per track, since AtGenfitter now writes them separately:
#   ke / theta / ex          back-extrapolated to the beam axis + CATIMA over the vertex gap
#   kefit / thetafit / exfit what the fit itself returned at the first measurement point
#
#   ./cache_dv135.sh [nparallel] [fitdir] [out]
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}

NPAR="${1:-4}"
STEM="${4-_multifit}"   # run-name stem: "_multifit" for the dv135 set, "" for gf_dt_matON/catima
FITDIR="${2:-/mnt/f/a1975/gf_dt_dv135/}"
OUT="${3:-/mnt/f/a1975/dt_kin_dv135.root}"
TMPD="$(dirname "$OUT")/.cache_$(basename "$OUT" .root)"
mkdir -p "$TMPD"
NUMS="0016 0017 0018 0019 0020 0021 0022 0023 0026 0027 0031 0032 0034 0036 0037 0038 0039 0040 0041 0042 0043 0044 0046 0048 0057 0058 0076 0077 0078 0079 0080 0082 0083 0084 0085 0086 0087 0088 0089 0091 0092 0095 0096 0097 0098 0102 0103"

one() {
  n="$1"; d="$2"; t="$3"; out="${t}/run_${n}.root"
  [ -s "$out" ] && { echo "[skip] run_$n"; return; }
  [ -s "${d}run_${n}${STEM}_genfitter_t.root" ] || { echo "[none] run_$n"; return; }
  root -b -l -q "ex_dt_a1975.C(\"run_${n}${STEM}\",\"${d}\",\"\",\"pid/triton_d2.json\",184.0,10.0,0.0,90.0,10.0,90.0,\"${out}\",\"/dev/null.png\")" \
      > "${t}/log_${n}.txt" 2>&1
  [ -s "$out" ] && echo "[done] run_$n" || echo "[FAIL] run_$n"
}
export -f one; export STEM
printf '%s\n' $NUMS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {} "$FITDIR" "$TMPD"
hadd -f "$OUT" "$TMPD"/run_*.root > "$TMPD/hadd.log" 2>&1
echo "=== cache: $OUT  ($(ls "$TMPD"/run_*.root 2>/dev/null | wc -l) runs) ==="
