#!/usr/bin/env bash
# Overnight: acceptance vs theta_cm for 14C(p,p'), ground state AND first excited state.
#
# Two samples, identical but for the residual excitation:
#   gs   Ex = 0      (elastic)
#   ex1  Ex = 6.094  (14C first excited state, 1-)
#
# FULL CM RANGE 2-178 deg on purpose. The production sims used 5-120 because beyond ~120 the
# recoil proton is too soft to reconstruct -- but for an ACCEPTANCE measurement that region is
# precisely what has to be measured (acceptance -> 0 there), so truncating the generator would
# build the answer into the input.
#
# The fit runs on the UNGATED reco: with the 3 cm beam hole the beam makes no hits at all, so
# the ungated reco is already essentially pure recoil protons, and keeping every event means the
# fit file has one entry per generated event -- which is what lets acceptance_C14.C pair truth
# with reconstruction by entry index. acceptance_C14.C refuses to run if the counts disagree.
set -eo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); SIM=$(cd "$HERE/.." && pwd)
REPO=$(cd "$SIM/../../../.." && pwd); UKF="$REPO/macro/Unpack_HDF5/a1954/UKF"
OUT="$HERE/acc"; mkdir -p "$OUT"
NEV=${1:-8000}
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
cd "$SIM"

run_one () {
  local TAG=$1 EX=$2
  echo "[$(date +%H:%M:%S)] ===== $TAG : residual Ex = $EX MeV ====="
  if [ ! -s "$OUT/${TAG}_sim.root" ]; then
    root -b -q -l "C14_pp_sim.C($NEV,2.0,178.0,\"TGeant4\",-28.5,\"$OUT/${TAG}_sim.root\",$EX)" > "$OUT/${TAG}_gen.log" 2>&1
    [ -s "$OUT/${TAG}_sim.root" ] || { echo "$TAG GEN_FAILED"; return 1; }
  fi
  echo "[$(date +%H:%M:%S)] $TAG reco (hole 30 mm, join mover)"
  if [ ! -s "$OUT/${TAG}_reco.root" ]; then
    root -b -q -l "run_reco_C14.C(\"$OUT/${TAG}_sim.root\",\"$OUT/${TAG}_reco.root\",\"ATTPC.a1954_C14_sim.par\",20,20,8,0,30.0,\"mover\")" > "$OUT/${TAG}_reco.log" 2>&1
    [ -s "$OUT/${TAG}_reco.root" ] || { echo "$TAG RECO_FAILED"; return 1; }
  fi
  echo "[$(date +%H:%M:%S)] $TAG UKF fit (bFieldSign -1, data convention)"
  root -b -q -l "$UKF/pipeline/fitUKF_C14.C(\"$TAG\",-1,\"proton\",-1,2.85,3.308e-5,\"\",\"$OUT/\",0.5,0.1,1,10,\"$OUT/\")" > "$OUT/${TAG}_fit.log" 2>&1
  [ -s "$OUT/${TAG}_ukf.root" ] || { echo "$TAG FIT_FAILED"; return 1; }
  echo "[$(date +%H:%M:%S)] $TAG acceptance"
  root -b -q -l "acceptance_C14.C(\"$OUT/${TAG}_sim.root\",\"$OUT/${TAG}_ukf.root\",\"$TAG\",$EX)" > "$OUT/${TAG}_acc.log" 2>&1
  grep -q "overall acceptance" "$OUT/${TAG}_acc.log" || { echo "$TAG ACC_FAILED"; return 1; }
  echo "[$(date +%H:%M:%S)] $TAG done"
}

run_one gs  0.0    || echo "GS STAGE FAILED"
run_one ex1 6.094  || echo "EX1 STAGE FAILED"

echo "===== SUMMARY ====="
for T in gs ex1; do
  echo "--- $T"; grep -A40 "overall acceptance" "$OUT/${T}_acc.log" 2>/dev/null | head -30
done
echo ACCEPTANCE_ALL_COMPLETED
