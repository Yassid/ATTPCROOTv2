#!/bin/bash
# UKF fit batch for a1954 12Be(p,p') — reads <run>_reco.root -> <run>_ukf.root
# Usage: ./fit_batch.sh "run_0142 run_0147 run_0150" [Nparallel] [particles] [bFieldSign]
# CPU-bound (reads the small _reco.root) -> 4-parallel is fine. Resumable.
RUNS="${1:-run_0142}"
NPAR="${2:-4}"
PARTS="${3:-proton}"
BSIGN="${4:--1}"
HERE="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="/home/yassid/a1954_Be12_reco"
LOGDIR="$OUTDIR/logs"; mkdir -p "$LOGDIR"

source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh
source /home/yassid/fair_install/ATTPCROOTv2/build/config.sh >/dev/null 2>&1

fit_one() {
  local run="$1"
  if [ ! -f "$OUTDIR/${run}_reco.root" ]; then echo "skip $run (no reco)"; return; fi
  if [ -f "$OUTDIR/${run}_ukf.root" ]; then echo "skip $run (ukf exists)"; return; fi
  echo "[$(date +%H:%M:%S)] ukf $run ..."
  root -b -q -l "$HERE/pipeline/fitUKF_Be12.C(\"$run\", -1, \"$PARTS\", $BSIGN)" > "$LOGDIR/${run}_ukf.log" 2>&1
  echo "[$(date +%H:%M:%S)] done $run -> exit $?"
}
export -f fit_one; export OUTDIR HERE LOGDIR PARTS BSIGN

echo "UKF fit runs: $RUNS  (${NPAR}-parallel, particles=$PARTS, bFieldSign=$BSIGN)"
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'fit_one "$@"' _ {}
echo "Fit batch complete."
