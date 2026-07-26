#!/bin/bash
# Chase the reco batch: as soon as a run is reported "done" in the reco batch log, slim it and
# run the full fit pass (gate -> UKF -> GENFIT) on it. Lets the whole 14C production finish
# unattended while reco is still chewing through the big runs.
#   ./follow_C14.sh "run_0065 run_0063 ..."
RUNS="${1:-}"
HERE="$(cd "$(dirname "$0")" && pwd)"
BATCHLOG=/home/yassid/a1954_C14_reco_hdb/batch.log
LOG=/home/yassid/a1954_C14_fit/follow.log
mkdir -p /home/yassid/a1954_C14_fit
for r in $RUNS; do
  # wait for reco to report this run done (and the file to stop growing)
  while ! grep -q "done $r " "$BATCHLOG" 2>/dev/null; do sleep 30; done
  [ -f "/home/yassid/a1954_C14_reco_hdb/${r}_reco.root" ] || { echo "$(date +%H:%M:%S) $r: no reco file" >> "$LOG"; continue; }
  [ -f "/home/yassid/a1954_C14_fit/${r}_genfit.root" ] && { echo "$(date +%H:%M:%S) $r: already fitted" >> "$LOG"; continue; }
  echo "$(date +%H:%M:%S) $r: slim+fit start" >> "$LOG"
  "$HERE/slim_cache_batch.sh" "$r" 1 >> "$LOG" 2>&1
  "$HERE/fitpipe_C14.sh" "$r" 2   >> "$LOG" 2>&1
  echo "$(date +%H:%M:%S) $r: done" >> "$LOG"
done
echo "$(date +%H:%M:%S) FOLLOW DONE" >> "$LOG"
