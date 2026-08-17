#!/usr/bin/env bash
# Build the (p,d) kinematics cache for the CATIMA material-effects production.
#
# Identical to how pd_kin.root was built -- cache_pd_run.C per run then hadd -- so the two
# caches differ ONLY by the fitter configuration and can be differenced directly. The one
# argument that changes is the file suffix: this production writes *_genfitter_pdcat.root
# while the existing one writes *_genfitter_pd.root.
#
# NOTHING IS OVERWRITTEN. Output is caches/pd_kin_catima.root, alongside the existing
# pd_kin.root, which is the input to the analysis kit and stays as it is.
#
# The IC gate is left at cache_pd_run's own default (tb 1000-1350) exactly as pd_kin.root was
# built, so the comparison is like-for-like. The ion-chamber value is STORED, not cut on, so
# the beam gate remains selectable at analysis time.
#
# WHY THE COMPLETION CHECK IS NOT `[ -s "$out" ]` (2026-08-17)
#   A reboot killed this script mid-write and left five per-run caches holding 276 bytes of a
#   300-byte ROOT header. Non-empty, so a restart called them "[have]" and skipped them, and the
#   hadd would have merged 79 runs while printing success -- a silently short cache is far worse
#   than a crashed job, because nothing downstream can tell. A run is now finished only when
#   pd/root_ok.C confirms the file opens, was closed cleanly, and holds a TTree; only then is a
#   .done marker written. Reruns cost one stat() per finished run, so the ROOT check is paid once.
#
#   Files predating the marker are validated once and adopted, so the good caches from an
#   interrupted run are kept rather than rebuilt.
#
#   The merge REFUSES to run if any run has a fit file but no valid cache, because that is the
#   exact condition that produces a short cache. FORCE=1 merges anyway, after printing what is
#   missing. Runs with no fit file at all are reported separately -- those are an expected skip,
#   not a failure.
#
#   ./cache_pd_catima.sh [nparallel]
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/UKF
# config.sh reads unset vars: source it BEFORE set -u or it kills the script silently
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
set -uo pipefail

NPAR="${1:-6}"
FORCE="${FORCE:-0}"
# Defaults are the BACK-EXTRAPOLATED production (_pdcatbx). The four paths move together, so the
# earlier no-back-extrapolation production is still rebuildable without editing this file:
#   GF=/mnt/f/a1975/reco_pd_catima/ SUF=_genfitter_pdcat \
#   WORK=/mnt/f/a1975/caches/.pdcat OUT=/mnt/f/a1975/caches/pd_kin_catima.root ./cache_pd_catima.sh
export GF="${GF:-/mnt/f/a1975/reco_pd_catima_bx/}"
export FRIB="${FRIB:-/mnt/f/a1975/reco/}"
export WORK="${WORK:-/mnt/f/a1975/caches/.pdcatbx}"
export OUT="${OUT:-/mnt/f/a1975/caches/pd_kin_catima_bx.root}"
export SUF="${SUF:-_genfitter_pdcatbx}"
mkdir -p "$WORK" /mnt/f/a1975/caches

# Echoes "VALID <entries>" or "INVALID <reason>"; empty output (a crashed ROOT) reads as invalid.
ok() {
  [ -s "$1" ] || { echo "INVALID missing"; return 0; }
  root -b -l -q "pd/root_ok.C(\"$1\")" 2>/dev/null | grep -E '^(VALID|INVALID) ' | tail -1
}
export -f ok

one() {
  r="$1"; run="run_${r}"
  out="${WORK}/${run}.root"

  # Cheap path: a marker means this run was already validated.
  [ -e "${out}.done" ] && [ -s "$out" ] && { echo "[have] $run"; return 0; }

  # A file with no marker is either pre-marker work or wreckage. Ask, do not assume.
  if [ -s "$out" ]; then
    v=$(ok "$out")
    case "$v" in
      VALID*) touch "${out}.done"; echo "[adopt] $run  (${v})"; return 0 ;;
      *)      echo "[stale] $run  (${v}) -- rebuilding"; rm -f "$out" "${out}.done" ;;
    esac
  fi

  [ -s "${GF}${run}${SUF}.root" ] || { echo "[nofit] $run"; return 0; }

  root -b -l -q "pp/cache_pd_run.C(\"${run}\",\"${out}\",\"$GF\",\"$FRIB\",\"$SUF\")" \
    > "${WORK}/log_${r}.txt" 2>&1

  v=$(ok "$out")
  case "$v" in
    "VALID 0") touch "${out}.done"; echo "[ok] $run  0 entries -- no tracks passed, check ${WORK}/log_${r}.txt" ;;
    VALID*)    touch "${out}.done"; echo "[ok] $run  ${v#VALID }" ;;
    *)         rm -f "$out"; echo "[FAIL] $run  (${v:-no answer}) see ${WORK}/log_${r}.txt" ;;
  esac
}
export -f one

RUNS=$(seq -f "%04g" 106 189)
echo "=== (p,d) CATIMA cache: $(echo $RUNS | wc -w) runs, $NPAR parallel ==="
printf '%s\n' $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}

# Account for every run before merging: a cache is only allowed to be short on purpose.
missing=""; noinput=""
for r in $RUNS; do
  if [ -e "${WORK}/run_${r}.root.done" ] && [ -s "${WORK}/run_${r}.root" ]; then continue; fi
  if [ -s "${GF}run_${r}${SUF}.root" ]; then missing="$missing $r"; else noinput="$noinput $r"; fi
done

n=$(ls "$WORK"/run_*.root 2>/dev/null | wc -l)
[ -n "$noinput" ] && echo "=== no fit file, expected skip:$noinput ==="
if [ -n "$missing" ]; then
  echo "=== INCOMPLETE: these runs have a fit file but no valid cache:$missing ==="
  if [ "$FORCE" != "1" ]; then
    echo "=== NOT MERGING. Fix them, or rerun with FORCE=1 to merge a knowingly short cache. ==="
    exit 1
  fi
  echo "=== FORCE=1: merging WITHOUT the runs above -- $OUT is short by $(echo $missing | wc -w) runs ==="
fi

echo "=== $n per-run caches; merging ==="
hadd -f "$OUT" "$WORK"/run_*.root > "$WORK/hadd.log" 2>&1 \
  && echo "=== cache: $OUT ===" || echo "=== HADD FAILED, see $WORK/hadd.log ==="
