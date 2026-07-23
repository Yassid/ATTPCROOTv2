#!/bin/bash
# Reco batch for a1954 12Be(p,p') — unpack + MultiFit PSA + clean + PRA -> <run>_reco.root
# Usage: ./reco_batch.sh "run_0142 run_0147 run_0150" [Nparallel]
# Resumable: skips runs whose _reco.root already exists.
# NOTE: raw HDF5 live on the external NSCL_e15250 drive; 2-parallel is the safe cap
#       (4-parallel thrashes the drive on the 5-9 GB files — see a1975 README).
RUNS="${1:-run_0142}"
NPAR="${2:-2}"
HERE="$(cd "$(dirname "$0")" && pwd)"
OUTDIR="/home/yassid/a1954_Be12_reco"
LOGDIR="$OUTDIR/logs"; mkdir -p "$LOGDIR"

source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh
source /home/yassid/fair_install/ATTPCROOTv2/build/config.sh >/dev/null 2>&1

run_one() {
  local run="$1"
  if [ -f "$OUTDIR/${run}_reco.root" ]; then echo "skip $run (exists)"; return; fi
  echo "[$(date +%H:%M:%S)] reco $run ..."
  root -b -q -l "$HERE/pipeline/unpackReco_Be12.C(\"$run\")" > "$LOGDIR/${run}_reco.log" 2>&1
  echo "[$(date +%H:%M:%S)] done $run -> exit $? ($(ls -la $OUTDIR/${run}_reco.root 2>/dev/null | awk '{printf "%.0f MB",$5/1e6}'))"
}
export -f run_one; export OUTDIR HERE LOGDIR

echo "Reco runs: $RUNS  (${NPAR}-parallel)"
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'run_one "$@"' _ {}
echo "Reco batch complete."
