#!/usr/bin/env bash
# Build the 16C(p,p) elastic kinematics cache over all 84 H2 runs.
#
# cache_pp_run.C's header describes being "driven over all runs by xargs -P4, then hadd into
# /tmp/pp_kin.root". That is how the cache behind the 316.4 mb^-1 elastic luminosity was made, and
# it is why that cache no longer exists: /tmp did not survive a reboot, and no script in the repo
# could rebuild it. An inline driver is not reproducible. This is.
#
# Output lands in /mnt/f/a1975/caches/ alongside the (p,d) and (d,t) caches, NOT in /tmp.
#
# Cost: one full pass over the FRIB files, ~88 GB, ~45 min at 7 shards. The ion chamber is not a
# stored variable -- see ANALYSIS_PROCEDURE.md section 3.
#
#   ./pp/cache_pp_all.sh [nparallel] [suffix]
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"; cd "$HERE"
NPAR="${1:-7}"
SUF="${2:-_genfitter_pphand}"
set +u; source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
WORK=/mnt/f/a1975/caches/.ppshards; mkdir -p "$WORK"; export WORK SUF HERE
one() {
  r="$1"; run="run_${r}"; out="$WORK/pp_${r}.root"
  [ -s "$out" ] && { echo "[have] $run"; return 0; }
  root -b -l -q "pp/cache_pp_run.C(\"$run\",\"$out\",\"/mnt/f/a1975/reco/\",\"$SUF\")" > "$WORK/log_${r}.txt" 2>&1
  if [ -s "$out" ]; then
    echo "[ok] $(grep -oE '^run_[0-9]+: cached.*' "$WORK/log_${r}.txt" | tail -1)"
    # the entry-count checks in cache_pp_run.C: surface ALL of them in the driver log, not just
    # the short-FRIB one. Greping for "FRIB SHORT" alone is how run_0157 and run_0158 stayed
    # invisible here for eight days -- their fit files were empty/truncated, which is a different
    # message, so the driver log looked clean.
    grep -aoE '(FRIB SHORT|FIT FILE EMPTY|more than the one)[^\\]*' "$WORK/log_${r}.txt" | sed "s/^/[BAD] $run: /"
  else
    echo "[FAIL] $run"
  fi
}
export -f one
echo "=== (p,p) elastic cache: runs 106-189, $NPAR shards, suffix $SUF ==="
seq -f "%04g" 106 189 | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
n=$(ls "$WORK"/pp_0*.root 2>/dev/null | wc -l); echo "=== $n run files; merging ==="
hadd -f /mnt/f/a1975/caches/pp_kin.root "$WORK"/pp_0*.root > "$WORK/hadd.log" 2>&1 \
  && echo "=== cache: /mnt/f/a1975/caches/pp_kin.root ===" || echo "=== HADD FAILED, see $WORK/hadd.log ==="
