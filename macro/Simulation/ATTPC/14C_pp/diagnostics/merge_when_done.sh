#!/usr/bin/env bash
# Waits for the 10 acceptance jobs, then merges each level across its 5 seeds.
set -eo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); SIM=$(cd "$HERE/.." && pwd)
REPO=$(cd "$SIM/../../../.." && pwd); OUT=/mnt/f/a1954_C14_acc
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
cd "$SIM"
while [ "$(ls $OUT/*.marker 2>/dev/null | wc -l)" -lt 10 ]; do
  pgrep -f "acc_parallel.sh|acc_batch.sh" >/dev/null || { echo "batch stopped with $(ls $OUT/*.marker 2>/dev/null | wc -l)/10 done"; break; }
  sleep 120
done
echo "=== markers: $(ls $OUT/*.marker 2>/dev/null | wc -l)/10"
root -b -q -l 'merge_acceptance.C("gs","1001,1002,1003,1004,1005",0.0)'   2>&1 | grep -vE "^Processing|^$" | tee "$OUT/merged_gs.txt"  | grep -E "MERGED|WARNING"
root -b -q -l 'merge_acceptance.C("ex1","1001,1002,1003,1004,1005",6.094)' 2>&1 | grep -vE "^Processing|^$" | tee "$OUT/merged_ex1.txt" | grep -E "MERGED|WARNING"
echo MERGE_COMPLETED
