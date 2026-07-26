#!/bin/bash
# Memory watchdog for the a1954 14C batches.
# Kills every ATTPCROOT root.exe job if NON-CACHE memory use crosses LIMIT% of the WSL cap.
# Page cache (which is what makes vmmem.exe look huge while streaming HDF5 off F:) is
# reclaimable and is NOT counted.
#   ./memwatch.sh [limit_percent] [poll_seconds]
LIM="${1:-80}"; POLL="${2:-30}"; LOG="/home/yassid/a1954_C14_reco_hdb/memwatch.log"
while true; do
  read -r TOT USED <<<"$(free -m | awk '/^Mem:/{print $2, $3}')"
  PCT=$((100 * USED / TOT))
  echo "$(date +%H:%M:%S) used=${USED}MB/${TOT}MB (${PCT}%)" >> "$LOG"
  if [ "$PCT" -ge "$LIM" ]; then
    echo "$(date +%H:%M:%S) OVER ${LIM}% -- killing reco/FRIB jobs" >> "$LOG"
    pkill -f "reco_hdb_C14.sh"; pkill -f "frib_C14_batch.sh"
    pkill -f "unpackReco_C14.C"; pkill -f "unpackFRIB_C14.C"
    exit 1
  fi
  sleep "$POLL"
done
