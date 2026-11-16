#!/usr/bin/env bash
# 2026-11-12: stop the local H2 reco if C: (holds the WSL ext4.vhdx) drops below MINFREE GB.
# Killing is safe: reco_tau_h2_local.sh skips runs with .done; partial runs are redone on relaunch.
MINFREE=${MINFREE:-20}; LOG=${GLOG:-/mnt/f/a1975/logs_h2_tctau/guard_cfree.log}; PAT=${GPAT:-reco_tctau_h2_local[.]sh}   # 2026-11-16: GPAT/GLOG for the D2 run
while pgrep -f "$PAT" >/dev/null; do
  free=$(df -BG --output=avail /mnt/d | tail -1 | tr -dc 0-9)
  echo "$(date '+%F %T') C: free ${free} GB" >> "$LOG"
  if [ "$free" -lt "$MINFREE" ]; then
    echo "$(date '+%F %T') BELOW ${MINFREE} GB -> killing H2 reco" >> "$LOG"
    pkill -f "$PAT"; pkill -f "unpackReco_multifit"; pkill -f "pid_width_dp"; pkill -x cp
    exit 1
  fi
  sleep 60
done
echo "$(date '+%F %T') reco finished, guard exits" >> "$LOG"
