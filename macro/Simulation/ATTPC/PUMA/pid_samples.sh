#!/bin/bash
# Prepare matched-momentum (p~375 MeV/c) labeled samples for K/pi dE/dx PID:
#   pi : branch-8 sim (existing attpcsim.root)
#   K  : branch-10 sim at E=0.620 GeV
# Both digitized + fit with species="both" (4 UKF hypotheses) so the fit-chi2 PID
# baseline is available alongside the dE/dx PID. Big digi outputs go to F:.
set -e
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC/PUMA
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh 2>/dev/null
SP=./data
export PUMA_OUT=/mnt/f/puma_sweep/
N=2000

echo "=== [$(date +%H:%M:%S)] preserve pi sim ==="
cp -f data/attpcsim.root data/attpcsim_pi.root

echo "=== [$(date +%H:%M:%S)] generate K sim (branch10, E=0.620, $N evt) ==="
root -b -q "PUMA_test10_sim.C($N, 0.620)" > $SP/pid_ksim.log 2>&1
mv -f data/attpcsim.root data/attpcsim_K.root
cp -f data/attpcsim_pi.root data/attpcsim.root   # restore default input

echo "=== [$(date +%H:%M:%S)] digi+fit(both) pi sample -> F: ==="
PUMA_IN=data/attpcsim_pi.root \
  root -b -q "run_digi_ukf_genfit_test8.C($N, 8.0, false, \"mfimpulse\", \"both\", 8)" \
  > $SP/pid_digi_pi.log 2>&1
mv -f /mnt/f/puma_sweep/output_digi_both8.root /mnt/f/puma_sweep/output_digi_pi_pid.root
echo "   pi -> $(ls -la /mnt/f/puma_sweep/output_digi_pi_pid.root | awk '{print $5}') bytes"

echo "=== [$(date +%H:%M:%S)] digi+fit(both) K sample -> F: ==="
PUMA_IN=data/attpcsim_K.root \
  root -b -q "run_digi_ukf_genfit_test8.C($N, 8.0, false, \"mfimpulse\", \"both\", 8)" \
  > $SP/pid_digi_K.log 2>&1
mv -f /mnt/f/puma_sweep/output_digi_both8.root /mnt/f/puma_sweep/output_digi_K_pid.root
echo "   K  -> $(ls -la /mnt/f/puma_sweep/output_digi_K_pid.root | awk '{print $5}') bytes"

echo "=== [$(date +%H:%M:%S)] samples ready ==="
