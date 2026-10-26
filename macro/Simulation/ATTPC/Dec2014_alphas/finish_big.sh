#!/bin/bash
# Digitise the already-generated 30000-event sim at the CALIBRATED threshold and add
# pile-up at the measured rate, giving the high-statistics efficiency sample.
#   thr 1000  -> simulated track density matches the data (scanned; 3 was 300x too low)
#   pile-up 0.10 -> the rate measured from the data (>=2 long tracks in a 12 deg cone)
cd /home/yassid/fair_install/ATTPCROOTv2_fr19port/macro/Simulation/ATTPC/Dec2014_alphas
source /home/yassid/fair_install/ATTPCROOTv2_fr19port/setup_fr19port.sh >/dev/null 2>&1
set -u
D=eff_big
echo "[$(date +%T)] digitising 30000 events at threshold 1000"
root -l -b -q "rundigi_sim_11p4.C(\"./$D/sim.root\",\"./$D/digi.root\",kFALSE,kFALSE,1000)" \
     > $D/digi.log 2>&1 || { echo "digi FAILED"; exit 1; }
du -sh $D/digi.root
echo "[$(date +%T)] merging"
root -l -b -q "merge_beam_reaction.C(\"$D/digi.root\",\"$D/merged.root\")" > $D/merge.log 2>&1
grep -E "merged events|hits/event" $D/merge.log
echo "[$(date +%T)] adding pile-up at the measured 10% rate"
root -l -b -q "add_pileup.C(\"$D/merged.root\",\"$D/digi.root\",\"$D/real.root\",0.10)" 2>&1 | grep -E "hits/event|events,"
echo "[$(date +%T)] BIG SAMPLE DONE"
