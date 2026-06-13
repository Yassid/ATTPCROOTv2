#!/bin/bash
# $1 = nIter. UKF deuteron per-cluster ref-track at nIter=$1 -> reco_pd/_rt$1
NIT=$1
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF
source build/config.sh >/dev/null 2>&1
cd macro/Unpack_HDF5/a1975/UKF
export RECO=/mnt/f/a1975/reco/ OUT=/mnt/f/a1975/reco_pd/ NIT
fit_one() {
  r=$1; run=run_$r; outf=${OUT}${run}_ukf_rt${NIT}.root
  if [ -f "$outf" ] && [ $(stat -c %s "$outf") -gt 1000000 ]; then echo ">>> $run done"; return; fi
  root -b -q "pipeline/fitUKF_a1975.C(\"$run\", -1, \"deuteron\", -1, 2.85, 9.0e-5, \"_rt${NIT}\", \"$RECO\", 0.5, 0.1, ${NIT}, 10, \"$OUT\", true, 4.0)" \
    > /tmp/ukf_rt${NIT}_${r}.log 2>&1
  echo "    done $run (rt${NIT})"
}
export -f fit_one
printf '%s\n' 0106 0107 0108 0109 0110 0111 0112 0113 | xargs -P2 -I{} bash -c 'fit_one {}'
echo "=== ALL DONE (ref-track nIter=${NIT}) ==="
