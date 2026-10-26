#!/bin/bash
# Run the whole chain in order, starting from the unpacked ROOT files.
#   ./run_all.sh
# Step 1 needs ROOT on the PATH; everything after that is pure Python.
set -e
cd "$(dirname "$0")"
python3 00_check.py
./01_dump.sh
python3 02_track.py
for s in 03_kinematics.py 04_selection.py 05_acceptance.py 06_cross_section.py; do
  echo; echo "=============================================================== $s"
  python3 "$s"
done
echo; echo "Done.  Plots are in out/"
