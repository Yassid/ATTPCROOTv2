#!/bin/bash
# Recover Round-1 param-scan ranking by running dres_eval.C over the leftover
# r1 fits (8-run subset 0106-0113). Results -> pd/dres_scan_results.txt (durable).
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF
source build/config.sh >/dev/null 2>&1
cd macro/Unpack_HDF5/a1975/UKF
RUNS="run_0106,run_0107,run_0108,run_0109,run_0110,run_0111,run_0112,run_0113"
OUT=pd/dres_scan_results.txt
: > "$OUT"
# label  suffix   (baseline first; variants change one knob each)
for pair in \
  "baseline_gd9 _d" \
  "gasDensity_7e-5 _r1_gd7" \
  "gasDensity_8e-5 _r1_gd8" \
  "gasDensity_10e-5 _r1_gd10" \
  "gasDensity_11e-5 _r1_gd11" \
  "measSigma_0.3 _r1_ms03" \
  "measSigma_0.75 _r1_ms075" \
  "momFrac_0.05 _r1_mf005" \
  "momFrac_0.2 _r1_mf02" \
  "nIter_2 _r1_ni2" ; do
  label=${pair%% *}; suf=${pair##* }
  echo ">>> $label ($suf)"
  root -b -q "pd/dres_eval.C(\"$RUNS\",\"$label\",\"$suf\")" 2>/dev/null | grep '^DRES' | tee -a "$OUT"
done
echo "=== DONE -> $OUT ==="
