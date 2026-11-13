#!/usr/bin/env bash
# Production for the levels_15C acceptance and response.
# One job per excitation energy; the macro's thetaCM is the RESIDUAL's CM angle, which is the same
# convention as thcmProj() in levels_15C/scripts/excm.h, so 2-80 deg covers the data's 10-70.
set -eo pipefail
# LC_ALL=C: a Spanish locale makes printf and bc expect DECIMAL COMMAS, so printf "%05.2f" 5.98
# fails with "invalid number" and the job tags come out wrong.
export LC_ALL=C
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
OUT=${OUT:-/home/yassid/Data/C15dd_sim}
NEV=${NEV:-20000}
NPAR=${NPAR:-10}
EBEAM=${EBEAM:-191.25}
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
mkdir -p "$OUT/logs"
# Ex grid: the six (p,p') resonances plus a spread for the efficiency curve
EXLIST=${EXLIST:-"0.00 2.00 3.00 4.43 4.93 5.98 6.68 7.38 8.28 10.00 12.00"}
echo "=== $(date +%H:%M:%S)  $NEV events per point, $NPAR parallel, Ebeam $EBEAM ==="
one(){
  local ex="$1"; local tag=$(printf "ex%05.2f" "$ex")
  local f="$OUT/${tag}_sim.root"
  [ -s "$f" ] && { echo "[$(date +%H:%M:%S)] $tag exists"; return 0; }
  local seed=$(( 9000 + $(printf '%.0f' "$(echo "$ex * 100" | bc -l)") ))
  local t0=$SECONDS
  root -b -q "$HERE/C15_dd_sim.C($NEV, 2., 80., \"TGeant4\", -28.5, \"$f\", $ex, $seed, $EBEAM)" \
       >"$OUT/logs/${tag}.log" 2>&1
  if [ -s "$f" ]; then echo "[$(date +%H:%M:%S)] $tag OK  $(du -h "$f"|cut -f1)  $(( (SECONDS-t0)/60 )) min"
  else echo "[$(date +%H:%M:%S)] $tag FAILED (see $OUT/logs/${tag}.log)"; fi
}
export -f one; export OUT NEV EBEAM HERE
printf '%s\n' $EXLIST | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "=== $(date +%H:%M:%S) done: $(ls -1 $OUT/ex*_sim.root 2>/dev/null|wc -l) files ==="
