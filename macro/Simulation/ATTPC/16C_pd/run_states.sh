#!/usr/bin/env bash
# Generate the 16C(p,d)15C simulation for the four states the analysis extracts.
#
#   gs   0.000  ground state, 1/2+
#   ex1  0.740  5/2+
#   ex2  3.100
#   ex3  4.660
#
# Transfer kinematics change with every one of them -- theta_lab(theta_cm) and the deuteron
# energy both move -- so an acceptance measured on one state does not describe another. That is
# why this runs per level rather than once.
#
# Generation only. Digitisation, reconstruction and fitting are separate stages; this is the
# cheap part (a few minutes per level) and the one that has to be right first.
#
#   ./run_states.sh [nEvents] [seed] [outDir]
set -eo pipefail
NEV=${1:-8000}
SEED=${2:-1001}
OUT=${3:-/mnt/f/a1975_C16_pd_sim}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../../../.." && pwd)
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
mkdir -p "$OUT"
cd "$HERE"

# tag:excitation
for CFG in "gs:0.0" "ex1:0.740" "ex2:3.100" "ex3:4.660"; do
  TAG=${CFG%%:*}; EX=${CFG##*:}
  J="${TAG}_s${SEED}"
  if [ -f "$OUT/$J.marker" ]; then
    echo "[$(date +%H:%M:%S)] $J already COMPLETED"
    continue
  fi
  echo "[$(date +%H:%M:%S)] $J: 15C at Ex = $EX MeV, $NEV events"
  root -b -q -l "C16_pd_sim.C($NEV,2.,178.,\"TGeant4\",-28.5,\"$OUT/${J}_sim.root\",$EX,$SEED)" \
       > "$OUT/${J}_gen.log" 2>&1
  # "the file exists" is not the same as "the job finished"; require both the seed to have been
  # applied and the macro to have reported success, or a killed job looks like a good one.
  [ -s "$OUT/${J}_sim.root" ] || { echo "  $J GEN_FAILED -- see $OUT/${J}_gen.log"; continue; }
  grep -q "RNG seed requested: $SEED" "$OUT/${J}_gen.log" || { echo "  $J SEED_NOT_APPLIED"; continue; }
  grep -q "Macro finished successfully" "$OUT/${J}_gen.log" || { echo "  $J DID_NOT_FINISH"; continue; }
  echo COMPLETED > "$OUT/$J.marker"
  echo "  $J done ($(du -h "$OUT/${J}_sim.root" | cut -f1))"
done
echo "[$(date +%H:%M:%S)] all four states generated in $OUT"
