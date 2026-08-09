#!/usr/bin/env bash
# GENFIT acceptance for an ARBITRARY residual excitation energy, one seed per invocation.
#
# acc_batch_genfit.sh cannot start a new level on its own: it requires the MC truth
# (<tag>_s<seed>_sim.root) to already exist in /mnt/f/a1954_C14_acc, because it was written to
# reuse the truth from the UKF acceptance pass so that the two fitters see identical events.
# For a level that has never been generated there is no such file and the job exits NO_TRUTH.
# This script fills that gap: it generates the truth if it is missing, then hands over.
#
# The reason it exists is the 8.5 and 9.4 MeV structures. Their yields were corrected with the
# 6.094 MeV acceptance, the only one available, and the gs-to-6.094 gradient says that
# extrapolation is worth about 2 percent. That is an argument, not a measurement, and this makes
# the measurement possible.
#
#   ./acc_level_genfit.sh <tag> <resEx> <seed> [nEvents]
#   ./acc_level_genfit.sh ex8 8.317 1001
#
# A marker is written only after the acceptance has actually printed a result, so its presence
# means COMPLETED rather than "a file appeared".
set -eo pipefail
TAG=$1; EX=$2; SEED=$3; NEV=${4:-8000}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); SIM=$(cd "$HERE/.." && pwd)
REPO=$(cd "$SIM/../../../.." && pwd)
SRC=/mnt/f/a1954_C14_acc            # where the shared MC truth lives
mkdir -p "$SRC"
J="${TAG}_s${SEED}"

set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
cd "$SIM"

if [ ! -s "$SRC/${J}_sim.root" ]; then
  echo "[$(date +%H:%M:%S)] $J: generating truth at Ex = $EX MeV, $NEV events"
  root -b -q -l "C14_pp_sim.C($NEV,2.0,178.0,\"TGeant4\",-28.5,\"$SRC/${J}_sim.root\",$EX,$SEED)" \
       > "$SRC/${J}_gen.log" 2>&1
  [ -s "$SRC/${J}_sim.root" ] || { echo "$J GEN_FAILED -- see $SRC/${J}_gen.log"; exit 1; }
  # the seed must actually have reached the generator, or every seed is the same sample
  grep -q "RNG seed requested: $SEED" "$SRC/${J}_gen.log" || { echo "$J SEED_NOT_APPLIED"; exit 1; }
  echo "[$(date +%H:%M:%S)] $J: truth done"
else
  echo "[$(date +%H:%M:%S)] $J: truth already present, reusing"
fi

exec "$HERE/acc_batch_genfit.sh" "$TAG" "$EX" "$SEED" "$NEV"
