#!/bin/bash
# genfit + CATIMA(H2) eloss on the 8-run subset, to isolate fitter vs eloss-model.
# Same as genfit_batch.sh but eloss=deuteron_H2_catima.txt -> reco_gf_catima/.
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF
source build/config.sh >/dev/null 2>&1
export GENFIT=/home/yassid/fair_install/GenFit
export LD_LIBRARY_PATH=$GENFIT/lib:$LD_LIBRARY_PATH
cd macro/Unpack_HDF5/a1975/UKF
mkdir -p /mnt/f/a1975/reco_gf_catima
export RECO=/mnt/f/a1975/reco/
export OUT=/mnt/f/a1975/reco_gf_catima/
fit_one() {
  r=$1; run=run_$r; outf=${OUT}${run}_genfit.root
  if [ -f "$outf" ] && [ $(stat -c %s "$outf") -gt 1000000 ]; then echo ">>> $run done, skip"; return; fi
  echo ">>> fitting $run (genfit+CATIMA) ..."
  root -b -q "pipeline/fitGenfit_a1975.C(\"$run\", -1, \"deuteron\", 2.85, 0.083147, \"$RECO\", \"\", \"$OUT\", \"deuteron_H2_catima.txt\")" \
    > /tmp/genfit_catima_${r}.log 2>&1
  echo "    done $run ($(grep -iE 'Real' /tmp/genfit_catima_${r}.log | tail -1))"
}
export -f fit_one
printf '%s\n' 0106 0107 0108 0109 0110 0111 0112 0113 | xargs -P2 -I{} bash -c 'fit_one {}'
echo "=== ALL DONE (genfit+CATIMA) ==="
