#!/usr/bin/env bash
# Generate the 16C(d,t)15C simulation for the five states the analysis extracts.
#
#   gs   0.000  ground state, 1/2+
#   ex1  0.740  5/2+
#   ex2  3.103
#   ex3  4.657
#   ex4  6.358  (the feature the data puts at 6.648 -- generated at the literature energy so the
#               acceptance is not conditioned on an assignment that is still open)
#
# Transfer kinematics change with every state, so an acceptance measured on one does not describe
# another. Ex also moves the theta_lab(theta_cm) mapping, and at Ex = 0 the tritons already sit at
# the 56 deg lab limit -- which is why the higher states are NOT a small correction here.
#
# CM RANGE 2-70 deg, and the number is a convention, not a truncation. AtTPC2Body ranges over the
# RESIDUAL's cm angle: asking for 110-178 returns ejectile angles of 5.7-70.7 deg. That is the same
# convention a1975_dt_excitation.py uses (theta_cm = pi - theta3_cm, "the residual recoils opposite
# the ejectile"), so 2-70 covers exactly the region the data occupies. Generating 2-178 instead
# would spend two thirds of the CPU on the high-KE branch, which has an ~840 mm cyclotron radius in
# a ~290 mm chamber and is not measurable.
#
# Generation only. Digitisation and reconstruction are separate stages.
#
#   ./run_states_dt.sh [nEvents] [seed] [outDir]
set -eo pipefail
NEV=${1:-20000}
SEED=${2:-3001}
OUT=${3:-/mnt/f/a1975_C16_dt_sim}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../../../.." && pwd)
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
mkdir -p "$OUT"
cd "$HERE"

for CFG in "gs:0.0" "ex1:0.740" "ex2:3.103" "ex3:4.657" "ex4:6.358"; do
  TAG=${CFG%%:*}; EX=${CFG##*:}
  J="${TAG}_s${SEED}"
  if [ -f "$OUT/$J.marker" ]; then
    echo "[$(date +%H:%M:%S)] $J already COMPLETED"
    continue
  fi
  echo "[$(date +%H:%M:%S)] $J: 15C at Ex = $EX MeV, $NEV events"
  root -b -q -l "C16_dt_sim.C($NEV,2.,70.,\"TGeant4\",-28.5,\"$OUT/${J}_sim.root\",$EX,$SEED)" \
       > "$OUT/${J}_gen.log" 2>&1
  # "the file exists" is not "the job finished": require the seed to have been applied AND the
  # macro to have reported success, or a killed job looks like a good one.
  [ -s "$OUT/${J}_sim.root" ] || { echo "  $J GEN_FAILED -- see $OUT/${J}_gen.log"; continue; }
  grep -q "RNG seed requested: $SEED" "$OUT/${J}_gen.log" || { echo "  $J SEED_NOT_APPLIED"; continue; }
  grep -q "Macro finished successfully" "$OUT/${J}_gen.log" || { echo "  $J DID_NOT_FINISH"; continue; }
  echo COMPLETED > "$OUT/$J.marker"
  echo "  $J done ($(du -h "$OUT/${J}_sim.root" | cut -f1))"
done
echo "[$(date +%H:%M:%S)] all five states generated in $OUT"
