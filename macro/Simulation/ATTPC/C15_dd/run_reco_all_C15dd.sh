#!/usr/bin/env bash
# Reconstruct every Ex point of the C15_dd production with the DATA chain at the VALIDATED gain.
# Gain 150000 (ATTPC.C15d_sim_g150000.par), not the data par's 10000: see README.md -- 10000 gives
# 56 hits/track against the data's 148 and silently deletes short tracks, which is what an
# acceptance measures.  holeR 20 mm, measured from the data's hit-radius ratio.
set -eo pipefail
export LC_ALL=C          # a Spanish locale breaks printf "%05.2f" and makes decimal-comma tags
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
OUT=${OUT:-/home/yassid/Data/C15dd_sim}
PAR=${PAR:-ATTPC.C15d_sim_g150000.par}
NPAR=${NPAR:-10}
NEV=${NEV:-0}            # 0 = all
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
mkdir -p "$OUT/logs"
one(){
  local f="$1"; local tag=$(basename "$f" _sim.root)
  local o="$OUT/${tag}_reco.root"
  [ -s "$o" ] && [ "$(stat -c%s "$o")" -gt 1000000 ] && { echo "[$(date +%H:%M:%S)] $tag exists"; return 0; }
  local t0=$SECONDS
  root -b -q "$HERE/run_reco_C15dd.C(\"$f\",\"$o\",\"$PAR\",20,15.0,7.5,$NEV,20.0)" \
       >"$OUT/logs/${tag}_reco.log" 2>&1 || true
  if [ -s "$o" ]; then echo "[$(date +%H:%M:%S)] $tag OK  $(du -h "$o"|cut -f1)  $(( (SECONDS-t0)/60 )) min"
  else echo "[$(date +%H:%M:%S)] $tag FAILED (see $OUT/logs/${tag}_reco.log)"; fi
}
export -f one; export OUT PAR NEV HERE
echo "=== $(date +%H:%M:%S)  reconstructing with $PAR, $NPAR parallel ==="
ls -1 "$OUT"/ex*_sim.root | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "=== $(date +%H:%M:%S) done: $(ls -1 $OUT/ex*_reco.root 2>/dev/null|wc -l) files ==="
