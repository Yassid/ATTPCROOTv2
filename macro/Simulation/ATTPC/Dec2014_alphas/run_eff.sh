#!/bin/bash
# Generate a large simulated sample for the tracking-efficiency measurement.
#   ./run_eff.sh 6000
cd "$(dirname "$0")"
N=${1:-6000}
source /home/yassid/fair_install/ATTPCROOTv2_fr19port/setup_fr19port.sh >/dev/null 2>&1
mkdir -p eff
echo "[$(date +%T)] simulating $N events"
root -l -b -q "He4He4_sim_el_11p4.C($N)" > eff/sim.log 2>&1 || { echo "sim failed"; exit 1; }
cp data11p4/attpcsim_in.root eff/sim.root
echo "[$(date +%T)] digitising"
root -l -b -q "rundigi_sim_11p4.C(\"./eff/sim.root\",\"./eff/digi.root\")" > eff/digi.log 2>&1 || { echo "digi failed"; exit 1; }
echo "[$(date +%T)] merging beam and reaction halves"
root -l -b -q "merge_beam_reaction.C(\"eff/digi.root\",\"eff/merged.root\")" > eff/merge.log 2>&1
grep -E "merged events|hits/event" eff/merge.log
echo "[$(date +%T)] done"
