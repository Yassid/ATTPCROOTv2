#!/bin/bash
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh 2>/dev/null
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
echo "===== REGEN START $(date) ====="
for r in 0300 0301 0302 0303 0304 0305; do
  echo "##### run_$r #####"
  root -l -b -q "unpackReco_multifit.C(\"run_$r\", 0, false)" 2>&1 | grep -iE "Done ->|Real time" | tail -2
done
echo "===== REGEN DONE $(date) ====="
