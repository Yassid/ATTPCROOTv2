#!/bin/bash
# Step 1 -- export every unpacked ROOT file to CSV.
#
# Needs ROOT on the PATH (any recent version; ATTPCROOT is NOT required -- the files
# carry their own StreamerInfo).  Check with:  root-config --version
#
#   ./01_dump.sh                 all runs
#   ./01_dump.sh 3               just the first 3, to try it out
cd "$(dirname "$0")"
command -v root >/dev/null || { echo "ROOT not on PATH. source thisroot.sh first."; exit 1; }
python3 - <<'PY' > /tmp/dec2014_paths.txt
import sys; sys.path.insert(0,'.')
import config, os
print(config.ROOT_DIR); print(config.CSV_DIR)
PY
ROOT_DIR=$(sed -n 1p /tmp/dec2014_paths.txt)
CSV_DIR=$(sed -n 2p /tmp/dec2014_paths.txt)
mkdir -p "$CSV_DIR"
N=${1:-0}
i=0
for f in "$ROOT_DIR"/alpha_run_*_hits.root; do
  [ -e "$f" ] || { echo "no ROOT files in $ROOT_DIR"; exit 1; }
  run=$(basename "$f" | sed 's/alpha_run_\([0-9]*\)_hits.root/\1/')
  csv="$CSV_DIR/hits_$run.csv"
  [ -s "$csv" ] && { echo "  skip $run (already dumped)"; continue; }
  root -l -b -q "lib/dump_hits.C(\"$f\",\"$csv\")" 2>&1 | grep -E "hits ->|ERROR"
  i=$((i+1)); [ "$N" -gt 0 ] && [ "$i" -ge "$N" ] && break
done
echo "  CSVs in $CSV_DIR"
