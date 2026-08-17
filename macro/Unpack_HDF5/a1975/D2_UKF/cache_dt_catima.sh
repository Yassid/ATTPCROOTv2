#!/usr/bin/env bash
# Build the (d,t) kinematics cache for the CATIMA production (prod_dt_catima.sh).
#
# Identical to cache_dt_dv1104.sh in EVERY selection -- same gate, same Ebeam, no cuts -- so the
# two caches differ only by the fitter configuration behind them and can be differenced directly.
# dt_kin_dv1104.root is material-effects OFF; this one is material effects ON with the CATIMA
# backend, the corrected dE/dx table, and the vertex-gap loss applied once instead of twice.
#
# The `ke` column comes from GetKinematicsXtr(), the back-extrapolated value -- which is exactly
# the slot the double count was inflating, so this cache is NOT comparable to any cache built
# before 2026-08-16 on a material-effects production.
#
# WHAT IS DIFFERENT FROM cache_dt_batch.sh, and why a separate script rather than an edit:
#
#   Ebeam   184.17 MeV, not 180.  11.5 MeV/u x 16.0147013 u.  Every cache built before
#           2026-08-14 used 180, and the cache's `ex` column is only valid at the energy it
#           was built with -- dt_flatness.C recomputes Ex from (ke, theta) precisely because
#           of that.  Getting this wrong looks exactly like a bad drift velocity.
#
#   inDir   /mnt/f/a1975/gf_dt_catima/, and the run string carries the _multifit tag because
#           ex_dt_a1975.C builds the filename as  inDir + run + "_genfitter_t" + suffix.
#           The IC file is still found: runTag is cut to the first 8 chars after "run_".
#
#   par     the geometry behind these fits is ATTPC.a1975_deuterium_dv1104.par, i.e.
#           dv 1.10424 cm/us, TBEntrance 560, ZPadPlane 971.7312 mm (NOT 1000 -- the drift
#           length follows from dv and the two TB anchors at 160 ns/TB).  z scale is
#           1.7668 mm/TB against 1.8182 before, so vertex-z cuts from the old production are
#           NOT transferable and must be re-derived.
#
# NO CUTS ARE APPLIED HERE beyond the PID gate (which genfit already applied anyway).
# chi2 -> 1e9, KE max -> 1e9, theta -> 0..180, and the IC amplitude is STORED rather than cut on.
# This is a NEW analysis: every selection is chosen on THIS data, not inherited. The old
# chi2 < 10 / KE < 90 / theta 10-90 defaults of ex_dt_a1975.C are exactly the kind of
# carried-over cut that quietly decides a result before anyone looks at it.
#
#   ./cache_dt_dv1104.sh [nparallel] [workdir]
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
NPAR="${1:-6}"
WORK="${2:-${DTWORK:-/mnt/f/a1975/caches/.cateloss}}"
OUT="${DTOUT:-/mnt/f/a1975/caches/dt_kin_cateloss.root}"
FIT="${DTFIT:-/mnt/f/a1975/gf_dt_cateloss/}"
EBEAM=184.17
mkdir -p "$WORK" /mnt/f/a1975/caches
# config.sh reads unset vars; under `set -u` it kills the script SILENTLY, zero-byte log,
# no error. Source it first, then turn -u on. This has now bitten twice.
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
set -uo pipefail
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}

NUMS="0016 0017 0018 0019 0020 0021 0022 0023 0026 0027 0031 0032 0034 0036 0037 0038 0039 0040 0041 0042 0043 0044 0046 0048 0057 0058 0076 0077 0078 0079 0080 0082 0083 0084 0085 0086 0087 0088 0089 0091 0092 0095 0096 0097 0098 0102 0103"

one() {
  n="$1"; work="$2"; fit="$3"; eb="$4"
  out="${work}/run_${n}.root"; log="${work}/log_${n}.txt"
  [ -s "$out" ] && { echo "[skip] run_$n"; return 0; }
  [ -s "${fit}run_${n}_multifit_genfitter_t.root" ] || { echo "[nofit] run_$n"; return 0; }
  root -b -l -q "ex_dt_a1975.C(\"run_${n}_multifit\",\"${fit}\",\"\",\"pid/triton_d2_dv1104.json\",${eb},1e9,0.0,1e9,0.0,180.0,\"${out}\",\"${work}/plot_${n}.png\")" \
     > "$log" 2>&1
  # the IC file is read over drvfs and fails INTERMITTENTLY under parallelism; a run that
  # silently lost its beam gate is worth ~1.7 MeV on Ex, so it is checked rather than assumed
  if grep -qi "no IC file" "$log"; then echo "[NO-IC] run_$n  <-- beam gate missing"; fi
  [ -s "$out" ] && echo "[done] run_$n $(date '+%H:%M:%S')" || echo "[FAIL] run_$n"
}
export -f one

printf '%s\n' $NUMS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {} "$WORK" "$FIT" "$EBEAM"

n=$(ls "$WORK"/run_*.root 2>/dev/null | wc -l)
echo "=== $n per-run caches; merging ==="
hadd -f "$OUT" "$WORK"/run_*.root > "$WORK/hadd.log" 2>&1 && echo "=== cache: $OUT ===" || echo "=== HADD FAILED, see $WORK/hadd.log ==="
grep -c "NO-IC" "$WORK"/../dv1104_cache.log 2>/dev/null || true
