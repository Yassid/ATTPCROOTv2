#!/bin/bash
# Scan UKF MomentumSigmaFrac vs 16C elastic FWHM, with MeasurementSigma fixed at 0.5.
source ~/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
cd ~/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/UKF/
OUT=/mnt/f/a1975/reco/
RUN=run_0106
NEV=${1:-15000}
MS=0.5
FRACS="0.1 0.2 0.3 0.5 0.7 1.0"

echo "=== UKF MomentumSigmaFrac scan ($RUN, $NEV evt, measSigma=$MS) ==="
for mf in $FRACS; do
  root -l -b -q "fitUKF_a1975.C(\"$RUN\", $NEV, \"proton\", -1, 2.85, 9.0e-5, \"_scan\", \"$OUT\", $MS, $mf)" \
    > /tmp/scan_mom_$mf.log 2>&1
  root -l -b -q "ex_eval.C(\"${OUT}${RUN}_ukf_scan.root\", \"${OUT}${RUN}_FRIB.root\", \"momFrac=$mf\")" \
    2>/dev/null | grep "^SCAN"
done
echo "=== SCAN DONE ==="
