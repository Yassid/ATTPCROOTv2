#!/usr/bin/env bash
# Build the 16C(d,p)17C kinematics cache from the CATIMA genfit production.
#
# Sibling of cache_dt_dv1104.sh, same structure, same `pk` tree, same per-run-then-hadd shape.
# Reads gf_dp_cateloss/run_NNNN_multifit_genfitter_p.root -> /mnt/f/a1975/caches/dp_kin_dv1104.root
#
# WHAT IS SPECIFIC TO THIS CHANNEL
#
#   EBEAM 184.25.  The (d,p) beam energy is carried THREE ways in this tree -- 184.25
#       (working_point.sh A1975_DT_EBEAM, MEASURED on these same 47 D2 runs and adopted for
#       (d,t)), 188.33 (C17_dp_fits README, 11.77 MeV/u, what the whole Spyral-era (d,p)
#       analysis used) and 192.0 (ex_dp_a1975.C's default, a purge-list value). 184.25 is taken
#       here because it is the only one measured on THIS data.
#       THIS IS NOT YET A CLOSED QUESTION, and it does not have to be closed before the cache is
#       built: the cache stores ke and theta per track, so Ex can be recomputed at any energy
#       afterwards without refitting -- which is exactly what dt_flatness.C does for (d,t). Only
#       the `ex`/`exfit` columns are energy-specific. Do not quote a 17C level off this cache
#       until the energy is settled; see the note in project_a1975_dp_analysis.
#
#   NO CUTS AT ALL.  chi2 -> 1e9, KE -> 1e9, theta -> 0..180, IC STORED not cut. The theta window
#       matters more here than in (d,t): the (d,p) transfer strength is BACKWARD, so the (d,t)
#       habit of 10-90 deg would throw away the physics. ex_dp_dv1104.C's own defaults are
#       already neutral, but they are passed explicitly so this file records the selection.
#
#   The gate is the LOOSE production gate. Tightening happens downstream off the stored
#       brho / dedx / sqrtdedx columns -- no refit needed, which is the whole point of caching
#       them (and is what (p,d) still cannot do, because cache_pd_run.C stores neither).
#
#   ./cache_dp_dv1104.sh [nparallel] [workdir]
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
NPAR="${1:-3}"   # 6 got the run OOM-killed on 2026-09-05; each job holds a ~400 MB fit file
WORK="${2:-/mnt/f/a1975/caches/.dp_dv1104}"
OUT="${OUT:-/mnt/f/a1975/caches/dp_kin_dv1104.root}"
FIT="${FIT:-/mnt/f/a1975/gf_dp_cateloss/}"
GATE="${GATE:-pid/proton_d2_dv1104.json}"
EBEAM="${EBEAM:-184.25}"
mkdir -p "$WORK" /mnt/f/a1975/caches
# config.sh reads unset vars; under `set -u` it kills the script SILENTLY, zero-byte log,
# no error. Source it first, then turn -u on.
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
set -uo pipefail
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}

[ -s "$GATE" ] || { echo "ERROR: PID gate '$GATE' missing. Draw it: root -l 'pid/run_gate_dp.C'"; exit 1; }
[ -d "$FIT" ]  || { echo "ERROR: no production at '$FIT'. Run ./prod_dp_catima.sh first."; exit 1; }

NUMS="0016 0017 0018 0019 0020 0021 0022 0023 0026 0027 0031 0032 0034 0036 0037 0038 0039 0040 0041 0042 0043 0044 0046 0048 0057 0058 0076 0077 0078 0079 0080 0082 0083 0084 0085 0086 0087 0088 0089 0091 0092 0095 0096 0097 0098 0102 0103"

one() {
  n="$1"; work="$2"; fit="$3"; eb="$4"; gate="$5"
  out="${work}/run_${n}.root"; log="${work}/log_${n}.txt"
  # SKIP GUARD -- must READ THE PRODUCT BACK, not test non-empty. This job was killed by the OOM
  # killer mid-write on 2026-09-05 and left 5 of 10 per-run caches as header-sized stubs; `[ -s ]`
  # calls those "present", so a naive restart skips them and the merged cache silently loses 5 runs
  # while reporting success. Exactly the failure pid/root_ok.C was written for.
  if [ -s "$out" ]; then
    vr=$(root -b -l -q "pid/root_ok.C(\"$out\")" 2>/dev/null | grep -E '^(VALID|INVALID)' | head -1)
    case "$vr" in
      VALID*) echo "[skip] run_$n  (${vr})"; return 0 ;;
      *)      echo "[redo] run_$n  (${vr:-unreadable}) -- removing stub"; rm -f "$out" ;;
    esac
  fi
  [ -s "${fit}run_${n}_multifit_genfitter_p.root" ] || { echo "[nofit] run_$n"; return 0; }
  root -b -l -q "ex_dp_dv1104.C(\"run_${n}_multifit\",\"${fit}\",\"\",\"${gate}\",${eb},1e9,0.0,1e9,0.0,180.0,\"${out}\",\"${work}/plot_${n}.png\")" \
     > "$log" 2>&1
  # the IC file is read over drvfs and fails INTERMITTENTLY under parallelism; a run that
  # silently lost its beam gate is worth ~1.7 MeV on Ex, so it is checked rather than assumed
  if grep -qi "NO IC after" "$log"; then echo "[NO-IC] run_$n  <-- beam gate missing"; fi
  [ -s "$out" ] && echo "[done] run_$n $(date '+%H:%M:%S')" || echo "[FAIL] run_$n"
}
export -f one

printf '%s\n' $NUMS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {} "$WORK" "$FIT" "$EBEAM" "$GATE"

n=$(ls "$WORK"/run_*.root 2>/dev/null | wc -l)
echo "=== $n per-run caches; merging ==="
hadd -f "$OUT" "$WORK"/run_*.root > "$WORK/hadd.log" 2>&1 && echo "=== cache: $OUT ===" || echo "=== HADD FAILED, see $WORK/hadd.log ==="
