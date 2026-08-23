#!/usr/bin/env bash
# Generate the 16C(p,t)14C simulation for the five 14C states the analysis would extract.
#
#   gs   0.000  ground state, 0+
#   ex1  6.094  1-
#   ex2  6.589  0+
#   ex3  6.728  3-
#   ex4  7.012  2+
#
# 14C has NOTHING between the ground state and 6.09, and the four excited states form a cluster
# this resolution (~0.6 MeV FWHM) cannot separate. They are generated separately anyway, because Ex
# moves the theta_lab(theta_cm) mapping and that lets the state dependence of the acceptance be
# MEASURED. (d,t) came out state-independent to 1%; whether that holds here is a result, not an
# assumption to inherit.
#
# CM RANGE 2-178, the FULL range, unlike run_states_dt.sh which truncates at 70. AtTPC2Body ranges
# over the RESIDUAL's cm angle; for (d,t) the backward branch has an ~840 mm cyclotron radius in a
# ~290 mm chamber and is unmeasurable, so truncating cost nothing. Here the radius peaks at about
# 298 mm near theta_cm 80-90 and falls again, so that branch is potentially IN acceptance -- and
# for the PID purity study it is the branch most worth having, being the stiff tracks most easily
# confused with something else.
#
# Generation only. Digitisation and reconstruction are separate stages.
#
#   ./run_states_pt.sh [nEvents] [seed] [outDir]
set -eo pipefail
NEV=${1:-20000}
SEED=${2:-4001}
OUT=${3:-/mnt/f/a1975_C16_pt_sim}
EBEAM=${EBEAM:-185.0}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../../../.." && pwd)
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
mkdir -p "$OUT" "$HERE/data"
cd "$HERE"

for CFG in "gs:0.0" "ex1:6.094" "ex2:6.589" "ex3:6.728" "ex4:7.012"; do
  TAG=${CFG%%:*}; EX=${CFG##*:}
  J="${TAG}_s${SEED}"
  if [ -f "$OUT/$J.marker" ]; then
    echo "[$(date +%H:%M:%S)] $J already COMPLETED"
    continue
  fi
  echo "[$(date +%H:%M:%S)] $J: 14C at Ex = $EX MeV, $NEV events, Ebeam $EBEAM"
  root -b -q -l "C16_pt_sim.C($NEV,2.,178.,\"TGeant4\",-28.5,\"$OUT/${J}_sim.root\",$EX,$SEED,$EBEAM)" \
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
