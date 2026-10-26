#!/bin/bash
# Steps 3-7 on the re-tracked data.  Run after run_all.sh completes.
cd /home/yassid/dec2014_package
export DEC2014_DATA=/home/yassid/dec2014_retrack
set -u
echo "[$(date +%H:%M:%S)] runs tracked: $(ls $DEC2014_DATA/tracks/*.pkl | wc -l)"
rm -f out/events_vd*.npz                 # the cache must not predate the re-tracking
python3 03_kinematics.py   2>&1 | tail -6
python3 04_selection.py    2>&1 | tail -8
python3 05_acceptance.py   2>&1 | tail -5
python3 06_cross_section.py 2>&1 | tail -5
python3 07_corrected_xsec.py 2>&1 | tail -6
echo "[$(date +%H:%M:%S)] done"
