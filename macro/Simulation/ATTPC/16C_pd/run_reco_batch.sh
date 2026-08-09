#!/usr/bin/env bash
# Digitise + reconstruct the generated (p,d) states. Slow: ~0.4 s/event, so ~50 min per state.
set -eo pipefail
SEED=${1:-s1001}; OUT=${2:-/mnt/f/a1975_C16_pd_sim}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); REPO=$(cd "$HERE/../../../.." && pwd)
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
cd "$HERE"
for TAG in gs ex1 ex2 ex3; do
  J="${TAG}_${SEED}"
  [ -f "$OUT/${J}_reco.marker" ] && { echo "$J reco already COMPLETED"; continue; }
  echo "[$(date +%H:%M:%S)] $J: reconstructing"
  root -b -q -l "run_reco_C16pd.C(\"$OUT/${J}_sim.root\",\"$OUT/${J}_reco.root\",\"ATTPC.a1975_C16_sim.par\",20,20,8,0,20.0,\"mover\")" \
       > "$OUT/${J}_reco.log" 2>&1
  # a written file is not a finished job; require the macro's own completion line
  if [ -s "$OUT/${J}_reco.root" ] && grep -q "sim reco done" "$OUT/${J}_reco.log"; then
    echo COMPLETED > "$OUT/${J}_reco.marker"; echo "  $J done"
  else
    echo "  $J RECO_FAILED -- see $OUT/${J}_reco.log"
  fi
done
echo "[$(date +%H:%M:%S)] all reco finished"
