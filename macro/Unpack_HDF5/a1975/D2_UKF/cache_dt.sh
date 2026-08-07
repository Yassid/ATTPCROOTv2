#!/usr/bin/env bash
# Build the (d,t) kinematics cache for one drift-velocity production.
#
# Ebeam is 184.0 MeV = 11.5 MeV/u, the correct beam energy. The older caches were built at 180
# and the free-Ebeam fits kept running to ~204, which is above the nominal incident energy and
# therefore an artefact (a momentum/B-field scale error being absorbed), not a beam energy.
#
# The cache carries BOTH energies per track, since AtGenfitter now writes them separately:
#   ke / theta / ex          back-extrapolated to the beam axis + CATIMA over the vertex gap
#   kefit / thetafit / exfit what the fit itself returned at the first measurement point
#
#   ./cache_dt.sh [nparallel] [fitdir] [out] [stem] [ebeam]
#     stem: "_multifit" for the dv sets (default), "" for gf_dt_matON / gf_dt_catima
#   e.g. ./cache_dt.sh 4 /mnt/f/a1975/gf_dt_dv1266/ /mnt/f/a1975/caches/dt_kin_dv1266.root
#
# RESUME SEMANTICS -- a marker means COMPLETED, never "a file exists".
# The old version resumed on `[ -s "$out" ]`, so a cache truncated by a crash, a killed job or
# a full disk was silently accepted and haddded into the final ntuple. A run is marked done
# only when ALL of these hold:
#   1. the macro printed its own completion line ("cached -> "), emitted after ntk->Write()
#   2. the log carries no segfault / ROOT stack trace
#   3. the output exists and clears a size floor
# The floor is deliberately TINY (1 kB, about the smallest a valid TNtuple file can be). An
# over-generous floor is the other half of this trap: a 100 MB reco floor once deleted the good
# output of run_0046/0103 because those runs only hold ~3k events. The completion LINE is the
# real signal; the floor only catches a zero-length or header-only file.
# A run that fails has its partial output DELETED, so the next pass genuinely redoes it.
set -o pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}

NPAR="${1:-4}"
FITDIR="${2:-/mnt/f/a1975/gf_dt_dv135/}"
OUT="${3:-/mnt/f/a1975/caches/dt_kin_dv135.root}"
STEM="${4-_multifit}"
EBEAM="${5:-184.0}"
TMPD="$(dirname "$OUT")/.cache_$(basename "$OUT" .root)"
mkdir -p "$TMPD"
NUMS="0016 0017 0018 0019 0020 0021 0022 0023 0026 0027 0031 0032 0034 0036 0037 0038 0039 0040 0041 0042 0043 0044 0046 0048 0057 0058 0076 0077 0078 0079 0080 0082 0083 0084 0085 0086 0087 0088 0089 0091 0092 0095 0096 0097 0098 0102 0103"

one() {
  n="$1"; d="$2"; t="$3"; out="${t}/run_${n}.root"; log="${t}/log_${n}.txt"
  [ -s "${out}.cachedone" ] && { echo "[skip] run_$n"; return 0; }
  [ -s "${d}run_${n}${STEM}_genfitter_t.root" ] || { echo "[none] run_$n"; return 0; }

  rm -f "$out"   # never let a previous partial masquerade as this pass's output
  root -b -l -q "ex_dt_a1975.C(\"run_${n}${STEM}\",\"${d}\",\"\",\"pid/triton_d2.json\",${EBEAM},10.0,0.0,90.0,10.0,90.0,\"${out}\",\"/dev/null.png\")" \
      > "$log" 2>&1

  if ! grep -q 'cached -> ' "$log"; then
    echo "[FAIL run_$n] macro never reached its cache write -- see $log"; rm -f "$out"; return 1
  fi
  if grep -qE 'segmentation violation|\*\*\* Break \*\*\*' "$log"; then
    echo "[FAIL run_$n] crashed after writing -- see $log"; rm -f "$out"; return 1
  fi
  if [ ! -s "$out" ] || [ "$(stat -c%s "$out" 2>/dev/null || echo 0)" -lt 1024 ]; then
    echo "[FAIL run_$n] output missing or header-only"; rm -f "$out"; return 1
  fi
  # candidate count is INFORMATIONAL, never a pass/fail test: a short run legitimately yields
  # few tritons, and gating on it would resurrect the size-floor trap in another form.
  ncand="$(grep -oE 'candidates [0-9]+' "$log" | tail -1 | awk '{print $2}')"
  touch "${out}.cachedone"
  echo "[done] run_$n  (${ncand:-?} candidates)"
}
export -f one; export STEM EBEAM
printf '%s\n' $NUMS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {} "$FITDIR" "$TMPD"

# --- loud accounting --------------------------------------------------------------------------
# A silent hadd over whatever happens to be on disk is how a partial production gets quoted as a
# complete one. Count what is present against what was asked for, and refuse to overwrite the
# final cache when runs are missing unless the caller explicitly forces it.
ndone=$(ls "$TMPD"/run_*.root.cachedone 2>/dev/null | wc -l)
nfit=0; for n in $NUMS; do [ -s "${FITDIR}run_${n}${STEM}_genfitter_t.root" ] && nfit=$((nfit+1)); done
nask=$(echo $NUMS | wc -w)
echo "=== runs asked $nask | fits present $nfit | cached OK $ndone ==="
if [ "$ndone" -lt "$nfit" ]; then
  echo "!!! $((nfit - ndone)) run(s) with a fit file FAILED to cache -- see $TMPD/log_*.txt"
  if [ "${FORCE:-0}" != "1" ]; then
    echo "!!! NOT writing $OUT (rerun to retry the failures, or FORCE=1 to hadd anyway)"; exit 1
  fi
  echo "!!! FORCE=1: hadding an INCOMPLETE set"
fi
if [ "$ndone" -eq 0 ]; then echo "!!! nothing cached, not writing $OUT"; exit 1; fi
hadd -f "$OUT" "$TMPD"/run_*.root > "$TMPD/hadd.log" 2>&1 || { echo "!!! hadd FAILED, see $TMPD/hadd.log"; exit 1; }
echo "=== cache: $OUT  ($ndone runs, Ebeam $EBEAM) ==="
