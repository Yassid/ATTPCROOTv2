#!/bin/bash
# UKF deuteron refit with reference-track warm start (nIter=2, refTrack=true),
# baseline params, eloss ON. -> reco_pd/_rt  for A/B vs baseline _d.
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF
source build/config.sh >/dev/null 2>&1
cd macro/Unpack_HDF5/a1975/UKF
export RECO=/mnt/f/a1975/reco/
export OUT=/mnt/f/a1975/reco_pd/
fit_one() {
  r=$1; run=run_$r; outf=${OUT}${run}_ukf_rt.root
  if [ -f "$outf" ] && [ $(stat -c %s "$outf") -gt 1000000 ]; then echo ">>> $run done, skip"; return; fi
  echo ">>> $run UKF ref-track ..."
  root -b -q "pipeline/fitUKF_a1975.C(\"$run\", -1, \"deuteron\", -1, 2.85, 9.0e-5, \"_rt\", \"$RECO\", 0.5, 0.1, 2, 10, \"$OUT\", true, 4.0)" \
    > /tmp/ukf_rt_${r}.log 2>&1
  echo "    done $run"
}
export -f fit_one
printf '%s\n' 0106 0107 0108 0109 0110 0111 0112 0113 | xargs -P2 -I{} bash -c 'fit_one {}'
echo "=== ALL DONE (UKF ref-track) ==="
