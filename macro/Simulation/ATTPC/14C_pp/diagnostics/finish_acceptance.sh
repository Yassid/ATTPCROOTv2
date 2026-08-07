#!/usr/bin/env bash
# Runs the acceptance step for both levels once run_acceptance.sh has produced the fits.
# Separate from the driver because acceptance_C14.C needs ROOT_INCLUDE_PATH (build/config.sh
# does not set it) and the driver was already running when that was discovered -- editing a
# running bash script in place can corrupt its execution, so it is left alone.
set -eo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); SIM=$(cd "$HERE/.." && pwd)
REPO=$(cd "$SIM/../../../.." && pwd); OUT="$HERE/acc"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
cd "$SIM"
for CFG in "gs:0.0" "ex1:6.094"; do
  T=${CFG%%:*}; EX=${CFG##*:}
  until [ -s "$OUT/${T}_ukf.root" ]; do
    pgrep -f run_acceptance.sh >/dev/null || { echo "driver gone before $T fit appeared"; exit 1; }
    sleep 60
  done
  sleep 20   # let the fit file close
  echo "[$(date +%H:%M:%S)] acceptance for $T (Ex = $EX)"
  root -b -q -l "acceptance_C14.C(\"$OUT/${T}_sim.root\",\"$OUT/${T}_ukf.root\",\"$T\",$EX)" 2>&1 \
    | grep -vE "^Processing|^$" | tee "$OUT/${T}_acceptance.txt" | grep -E "overall|ENTRY|no truth"
done
echo ACCEPTANCE_REPORT_COMPLETED
