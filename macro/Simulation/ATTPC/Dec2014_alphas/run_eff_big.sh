#!/bin/bash
# High-statistics sample for the tracking efficiency.
#
# The efficiency correction is limited by how many simulated reactions land in the
# 90-degree band: 3000 reactions gave only ~484 there, which is why the migration-aware
# bin-by-bin correction was unusable.  This runs 5x that.
#
# Both bulk branches are OFF (drifted electrons ~90% of the volume, pad traces most of the
# rest).  AtPSAtask consumes them in memory within the same run, so nothing downstream
# needs them, and digi.root drops from 50 GB to a few tens of MB.
#   ./run_eff_big.sh 30000
cd "$(dirname "$0")"
source /home/yassid/fair_install/ATTPCROOTv2_fr19port/setup_fr19port.sh >/dev/null 2>&1
set -u
N=${1:-30000}
D=eff_big
mkdir -p $D
echo "[$(date +%T)] simulating $N events"
root -l -b -q "He4He4_sim_el_11p4.C($N)" > $D/sim.log 2>&1 || { echo "sim FAILED"; exit 1; }
cp data11p4/attpcsim_in.root $D/sim.root
echo "[$(date +%T)] digitising (electrons and traces NOT persisted)"
root -l -b -q "rundigi_sim_11p4.C(\"./$D/sim.root\",\"./$D/digi.root\",kFALSE,kFALSE)" \
     > $D/digi.log 2>&1 || { echo "digi FAILED"; exit 1; }
du -sh $D/digi.root
echo "[$(date +%T)] merging beam and reaction halves"
root -l -b -q "merge_beam_reaction.C(\"$D/digi.root\",\"$D/merged.root\")" > $D/merge.log 2>&1
grep -E "merged events|hits/event" $D/merge.log
echo "[$(date +%T)] ALL DONE"
