#!/bin/bash
# Momentum scan for PUMA pi+/pi- (max PSA + ring clustering). For each beam energy
# it sims -> digi+UKF/GENFIT -> compares to truth, collecting momentum/angle/charge
# resolution vs |p|. Writes /tmp/momscan.txt (p_MeV  ukf_pbias ukf_psig ukf_theta
# ukf_charge  gf_pbias gf_psig gf_theta gf_charge).
source ~/fair_install/ATTPCROOTv2/build/config.sh 2>/dev/null
cd ~/fair_install/ATTPCROOTv2/macro/Simulation/ATTPC/PUMA
N=${1:-300}
OUT=/tmp/momscan.txt
: > $OUT
# (p_GeV  E_GeV) pairs; E = sqrt(p^2 + m_pi^2), m_pi=0.13957
for PE in "0.150 0.2046" "0.250 0.2863" "0.375 0.4001" "0.500 0.5191" "0.700 0.7138" "0.900 0.9107"; do
  P=$(echo $PE | cut -d' ' -f1); E=$(echo $PE | cut -d' ' -f2)
  echo ">>> momentum scan point p=$P GeV (E=$E GeV), N=$N"
  root -b -q "PUMA_test8_sim.C($N,$E)" > /tmp/ms_sim.log 2>&1
  root -b -q "run_digi_ukf_genfit_test8.C($N, 8.0, false, \"max\", \"pi\", 8, 16, 256, 0.0, true, true, false, 20)" > /tmp/ms_digi.log 2>&1
  mv -f data/output_digi_both8.root data/scan_p${P}.root
  # compare emits two blocks (UKF then GENFIT); grab p-res + theta + charge from each
  CMP=$(root -b -q "compare_ukf_genfit_test8.C(\"pi\", $E, \"./data/scan_p${P}.root\", \"./data/attpcsim.root\")" 2>/dev/null)
  # parse: lines "p resolution : median +X% sigma_IQR Y%", "theta ... median M sigma_IQR S", "charge-sign accuracy: Z%"
  readarray -t PBIAS < <(echo "$CMP" | grep -oE "median [+-]?[0-9.]+%  sigma_IQR" | grep -oE "[+-]?[0-9.]+" )
  readarray -t PSIG  < <(echo "$CMP" | grep -oE "sigma_IQR [0-9.]+%" | grep -oE "[0-9.]+")
  readarray -t THSIG < <(echo "$CMP" | grep -A1 "theta" | grep -oE "sigma_IQR [0-9.]+" | grep -oE "[0-9.]+")
  readarray -t CHG   < <(echo "$CMP" | grep -oE "charge-sign accuracy: [0-9.]+%" | grep -oE "[0-9.]+")
  PMEV=$(echo "$P*1000" | bc)
  echo "$PMEV ${PBIAS[0]:-0} ${PSIG[0]:-0} ${THSIG[0]:-0} ${CHG[0]:-0} ${PBIAS[1]:-0} ${PSIG[1]:-0} ${THSIG[1]:-0} ${CHG[1]:-0}" >> $OUT
  echo "   -> $(tail -1 $OUT)"
done
echo "MOMSCAN_DONE"; cat $OUT
