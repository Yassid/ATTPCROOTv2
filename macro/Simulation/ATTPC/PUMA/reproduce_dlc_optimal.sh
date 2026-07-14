#!/bin/bash
# reproduce_dlc_optimal.sh — reproduce the OPTIMAL PUMA reconstruction with DLC on.
#
#   Winning chain: AtPSAMultiFit(impulse) -> HDBSCAN -> per-ring charge-weighted
#   clustering -> UKF & GENFIT, DLC 1.35 MOhm/sq, Cu/Al material correction.
#
#   Expected @375 MeV/c (pi+pi-):
#       GENFIT :  bias ~0% ,  sigma ~7%   (2x better than no-DLC; sub-pad gain realised)
#       UKF    :  bias ~-6%,  sigma ~15%  (recovers no-DLC; residual bias calibratable)
#
#   Usage:  bash reproduce_dlc_optimal.sh [nEvents] [E_GeV]
#           nEvents default 500 ; E_GeV default 0.4001 (= |p|=375 MeV/c pion)
set -e
source ~/fair_install/ATTPCROOTv2/build/config.sh 2>/dev/null
cd ~/fair_install/ATTPCROOTv2/macro/Simulation/ATTPC/PUMA
N=${1:-500}
E=${2:-0.4001}

echo "[1/3] simulate $N back-to-back pi+pi-  (E=$E GeV/pion) ..."
root -b -q "PUMA_test8_sim.C($N,$E)" > /tmp/rdo_sim.log 2>&1
cp -f data/attpcsim.root data/opt_sim.root

echo "[2/3] optimal digi + reco + fit (DLC on, MultiFit, HDBSCAN, per-ring, GENFIT) ..."
root -b -q "reco_dlc_optimal.C($N,\"data/opt_sim.root\",\"data/opt_reco.root\")" > /tmp/rdo_reco.log 2>&1

echo "[3/3] UKF vs GENFIT resolution:"
root -b -q "compare_ukf_genfit_test8.C(\"pi\",$E,\"./data/opt_reco.root\",\"./data/opt_sim.root\")" 2>&1 \
  | sed 's/\x1b\[[0-9;]*m//g' | grep -E "===|p resolution|theta|charge|vertex dz" | grep -v Recon
echo ""
echo "DONE. Winning chain reproduced (GENFIT should be ~unbiased, sigma ~7%)."
