#!/bin/bash
# Fit the 8-run subset (0106-0113) with modified GenFit on the SAME _reco.root the
# UKF used, so dres_eval can compare fitter resolution head-to-head. Resumable
# (skips runs whose _genfit.root already exists >1MB). Output -> reco_gf/.
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF
source build/config.sh >/dev/null 2>&1
export GENFIT=/home/yassid/fair_install/GenFit
export LD_LIBRARY_PATH=$GENFIT/lib:$LD_LIBRARY_PATH
cd macro/Unpack_HDF5/a1975/UKF
mkdir -p /mnt/f/a1975/reco_gf
export RECO=/mnt/f/a1975/reco/
export OUT=/mnt/f/a1975/reco_gf/
fit_one() {
  r=$1; run=run_$r; outf=${OUT}${run}_genfit.root
  if [ -f "$outf" ] && [ $(stat -c %s "$outf") -gt 1000000 ]; then
    echo ">>> $run already done, skip"; return
  fi
  echo ">>> fitting $run with genfit ..."
  root -b -q "pipeline/fitGenfit_a1975.C(\"$run\", -1, \"deuteron\", 2.85, 0.083147, \"$RECO\", \"\", \"$OUT\")" \
    > /tmp/genfit_${r}.log 2>&1
  echo "    done $run -> $outf ($(grep -iE 'Real' /tmp/genfit_${r}.log | tail -1))"
}
export -f fit_one
# 4-parallel: genfit fits are CPU-bound and read the smaller _reco.root (page-cached),
# so they parallelize cleanly (unlike the HDF5 unpack which is drvfs-I/O-capped at 2).
# 4 x ~0.8 GB RSS = ~3 GB, well under the memory cap on this 8-core box.
printf '%s\n' 0106 0107 0108 0109 0110 0111 0112 0113 | xargs -P4 -I{} bash -c 'fit_one {}'
echo "=== ALL DONE ==="
