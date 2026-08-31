#!/usr/bin/env bash
# Wait for the 300 torr Magboltz B scan, build one par per field from it, then run the g.s.
# field comparison. Single core end to end.
set -eo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCAN=/home/yassid/attpc_dv_3He/out300
LOG=/mnt/f/ar46_3hed_bcompare.log
say() { echo "[$(date +%F' '%H:%M:%S)] CHAIN: $*" | tee -a "$LOG"; }

say "waiting for the Magboltz scan at 140 V/cm (v_d 1.32, 512 tb @ 6 MHz)"
until grep -q "op140 scan complete" "$SCAN/op140scan.log" 2>/dev/null; do
  if ! pgrep -f "[r]un_300torr_op140" >/dev/null && ! grep -q "op140 scan complete" "$SCAN/op140scan.log" 2>/dev/null; then
    say "SCAN DIED without completing -- stopping"; exit 1
  fi
  sleep 60
done
say "scan done; building pars"
"$DIR/make_3Hed_pars.sh" 2>&1 | tee -a "$LOG"
say "launching the field comparison"
exec "$DIR/run_3Hed_bcompare.sh"
