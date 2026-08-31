#!/bin/bash
# Refit every 46Ar(3He,d) sample with the post-2026-08-29 AtGenfitter: vertex-end measurement
# order, measSigma 0.6, matEffects + CATIMA, backward seed fix, bField NEGATIVE.
#
# SINGLE CORE ON PURPOSE -- the 14C(d,p) campaign owns the rest of the machine.
# Only the AT-TPC pad plane is refitted: the 2 mm configurations are parked because the beam hole
# has to stay for a heavy-fragment telescope.
#
# Resumable on EVIDENCE the fit finished ("Done." in the log AND a non-empty output), not on the
# file existing -- a killed root leaves a valid-looking truncated file behind.
set -eo pipefail
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
cd "$REPO"
export VMCWORKDIR=$PWD
source build/config.sh >/dev/null 2>&1

TAGS="gs_s3001 gs_s3002 360_s3011 360_s3012 2020_s3021 2020_s3022"
for cfg in "/mnt/f/ar46_3hed/:-2.85" "/mnt/f/ar46_3hed_B38/:-3.80"; do
  DIR=${cfg%%:*}; BF=${cfg##*:}
  for t in $TAGS; do
    OUT="$DIR${t}_genfitter_d.root"; LOG="$DIR${t}_fit.log"
    if grep -q "^.\[1;32mDone" "$LOG" 2>/dev/null && [ -s "$OUT" ]; then
      echo "[$(date +%H:%M:%S)] $DIR$t already fitted, skipping"; continue
    fi
    echo "[$(date +%H:%M:%S)] fitting $t  B=$BF  ($DIR)"
    root -b -q -l "macro/Simulation/ATTPC/46Ar_3Hed/fitGenfitter_Ar46.C(\"$t\",-1,\"$DIR\",\"\",\"\",$BF)" > "$LOG" 2>&1
    grep -q "Done" "$LOG" || { echo "$t FIT_FAILED"; exit 1; }
    [ -s "$OUT" ] || { echo "$t NO_OUTPUT"; exit 1; }
    # assert the material model really engaged: an inert arm must not pass for a valid one
    grep -q "dE/dx from CATIMA" "$LOG" || { echo "$t CATIMA_NOT_ENABLED"; exit 1; }
  done
done
echo "[$(date +%H:%M:%S)] refit complete"
