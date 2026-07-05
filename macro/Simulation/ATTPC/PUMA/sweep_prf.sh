#!/bin/bash
# PRF charge-dispersion sweep for PUMA branch-8 pions, grounded in SPSC-P-361 sec 7.4.1.
# Telegraph dispersion sigma = sqrt(2t/RC); pad pitch dl~2.3mm azimuthal -> sweep sigma
# from single-pad (0) to ~3-pad spread (1.5mm). Ring-centroiding ON for sigma>0.
# Outputs (~GB each) go to F: (PUMA_OUT) to avoid growing the WSL vhdx on C:.
set -u
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC/PUMA
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh 2>/dev/null
SP=/tmp/claude-1000/-home-yassid/18bdf430-6b82-49d5-b6bd-2bbc0fa19eb0/scratchpad
SWDIR=/mnt/f/puma_sweep
mkdir -p "$SWDIR"
export PUMA_OUT="$SWDIR/"
N=${1:-5000}

run_point () {
  local tag=$1 sig=$2 ring=$3
  local out="$SWDIR/output_digi_${tag}.root"
  echo "=== [$(date +%H:%M:%S)] point $tag : prfSigma=$sig ring=$ring  N=$N ==="
  root -b -q "run_digi_ukf_genfit_test8.C($N, 8.0, false, \"mfimpulse\", \"pi\", 8, 16, 256, $sig, $ring)" \
       > $SP/sweep_${tag}.log 2>&1
  mv -f "$SWDIR/output_digi_both8.root" "$out" 2>/dev/null
  ls -la "$out" 2>/dev/null | awk '{print "   -> "$5" bytes  "$9}'
}

# remaining physical dispersion points (base/ring0/0.3mm already on F:)
run_point sw_r05    0.5 true
run_point sw_r08    0.8 true
run_point sw_r12    1.2 true
run_point sw_r15    1.5 true
echo "=== [$(date +%H:%M:%S)] sweep complete ==="
