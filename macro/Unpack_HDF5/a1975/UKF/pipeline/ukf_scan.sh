#!/bin/bash
# Scan UKF MeasurementSigma vs 16C elastic-peak resolution (FWHM).
# Refits run_0106 (subset) for each value, suffix _scan so production is untouched.
source ~/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
cd ~/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/UKF/
OUT=/mnt/f/a1975/reco/
RUN=run_0106
NEV=${1:-8000}
shift 2>/dev/null
SIGMAS="${@:-0.5 0.75 1.0 1.5 2.0 3.0}"

echo "=== UKF MeasurementSigma scan ($RUN, $NEV events) ==="
for ms in $SIGMAS; do
  root -l -b -q "fitUKF_a1975.C(\"$RUN\", $NEV, \"proton\", -1, 2.85, 9.0e-5, \"_scan\", \"$OUT\", $ms)" \
    > /tmp/scan_fit_$ms.log 2>&1
  root -l -b -q "ex_eval.C(\"${OUT}${RUN}_ukf_scan.root\", \"${OUT}${RUN}_FRIB.root\", \"measSigma=$ms\")" \
    2>/dev/null | grep "^SCAN"
done
echo "=== SCAN DONE ==="
