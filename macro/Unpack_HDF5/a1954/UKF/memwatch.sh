#!/bin/bash
# Memory AND host-disk watchdog for the a1954 14C batches.
#
# Memory: kills every ATTPCROOT root.exe job if NON-CACHE memory use crosses LIMIT% of the
# WSL cap. Page cache (which is what makes vmmem.exe look huge while streaming HDF5 off F:)
# is reclaimable and is NOT counted.
#
# Disk: watches the WINDOWS C: drive, not `/`. The ext4 filesystem lives in a sparse
# ext4.vhdx on C:, so `df /` reports the guest view and stays at hundreds of GB free even
# as the host fills. On 2026-07-26 C: hit 100% and WSL was killed outright at 02:35:36,
# truncating three in-flight recos (run_0055/0058/0060). Guest df showed 814 GB "free"
# the whole time. Kill the jobs before that happens again.
#   ./memwatch.sh [limit_percent] [poll_seconds] [min_C_free_GB]
LIM="${1:-80}"; POLL="${2:-30}"; MINC="${3:-3}"
LOG="/home/yassid/a1954_C14_reco_hdb/memwatch.log"
kill_jobs(){
  pkill -f "reco_hdb_C14.sh"; pkill -f "frib_C14_batch.sh"; pkill -f "follow_C14.sh"
  pkill -f "unpackReco_C14.C"; pkill -f "unpackFRIB_C14.C"
}
while true; do
  read -r TOT USED <<<"$(free -m | awk '/^Mem:/{print $2, $3}')"
  PCT=$((100 * USED / TOT))
  CFREE=$(df -BG --output=avail /mnt/c 2>/dev/null | tail -1 | tr -dc '0-9')
  CFREE="${CFREE:-999}"
  echo "$(date +%H:%M:%S) used=${USED}MB/${TOT}MB (${PCT}%) C_free=${CFREE}G" >> "$LOG"
  if [ "$PCT" -ge "$LIM" ]; then
    echo "$(date +%H:%M:%S) OVER ${LIM}% MEM -- killing reco/FRIB jobs" >> "$LOG"
    kill_jobs; exit 1
  fi
  if [ "$CFREE" -le "$MINC" ]; then
    echo "$(date +%H:%M:%S) C: DOWN TO ${CFREE}G (<= ${MINC}G) -- killing reco/FRIB jobs" >> "$LOG"
    kill_jobs; exit 2
  fi
  sleep "$POLL"
done
