#!/usr/bin/env bash
# One acceptance job: generate -> reco -> fit -> acceptance, for a given level and RNG SEED.
# Output goes to the external drive F: (/mnt/f/a1954_C14_acc), because C: is down to ~70 GB and
# each job writes ~2.7 GB. Called by acc_parallel.sh; runnable standalone.
#
#   ./acc_batch.sh <tag> <resEx> <seed> [nEvents]
#
# The marker is written ONLY after acceptance_C14.C has actually printed a result, so
# "marker exists" means COMPLETED -- never merely "a file appeared".
set -eo pipefail
TAG=$1; EX=$2; SEED=$3; NEV=${4:-8000}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); SIM=$(cd "$HERE/.." && pwd)
REPO=$(cd "$SIM/../../../.." && pwd); UKF="$REPO/macro/Unpack_HDF5/a1954/UKF"
OUT=${ACC_OUT:-/mnt/f/a1954_C14_acc}
mkdir -p "$OUT"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
cd "$SIM"
J="${TAG}_s${SEED}"
[ -f "$OUT/$J.marker" ] && { echo "$J already COMPLETED"; exit 0; }

# Reuse a complete generation if one is already there -- regenerating costs ~1 min/job and the
# samples were already verified to differ per seed.
if [ ! -s "$OUT/${J}_sim.root" ]; then
  root -b -q -l "C14_pp_sim.C($NEV,2.0,178.0,\"TGeant4\",-28.5,\"$OUT/${J}_sim.root\",$EX,$SEED)" > "$OUT/${J}_gen.log" 2>&1
  [ -s "$OUT/${J}_sim.root" ] || { echo "$J GEN_FAILED"; exit 1; }
  grep -q "RNG seed requested: $SEED" "$OUT/${J}_gen.log" || { echo "$J SEED_NOT_APPLIED"; exit 1; }
else
  echo "$J: reusing existing generation"
fi

root -b -q -l "run_reco_C14.C(\"$OUT/${J}_sim.root\",\"$OUT/${J}_reco.root\",\"ATTPC.a1954_C14_sim.par\",20,20,8,0,30.0,\"mover\")" > "$OUT/${J}_reco.log" 2>&1
[ -s "$OUT/${J}_reco.root" ] || { echo "$J RECO_FAILED"; exit 1; }

root -b -q -l "$UKF/pipeline/fitUKF_C14.C(\"$J\",-1,\"proton\",-1,2.85,3.308e-5,\"\",\"$OUT/\",0.5,0.1,1,10,\"$OUT/\")" > "$OUT/${J}_fit.log" 2>&1
[ -s "$OUT/${J}_ukf.root" ] || { echo "$J FIT_FAILED"; exit 1; }

root -b -q -l "acceptance_C14.C(\"$OUT/${J}_sim.root\",\"$OUT/${J}_ukf.root\",\"$J\",$EX)" > "$OUT/${J}_acc.log" 2>&1
grep -q "overall acceptance" "$OUT/${J}_acc.log" || { echo "$J ACC_FAILED"; exit 1; }
mv -f "$HERE/acceptance_${J}.root" "$OUT/" 2>/dev/null || true
rm -f "$HERE/acceptance_${J}.png"

# the reco is the bulk of the space and nothing downstream needs it once the fit exists
rm -f "$OUT/${J}_reco.root"
echo COMPLETED > "$OUT/$J.marker"
echo "$J done: $(grep 'overall acceptance' "$OUT/${J}_acc.log")"
