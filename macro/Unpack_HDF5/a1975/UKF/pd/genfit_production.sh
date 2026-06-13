#!/bin/bash
# GenFit (p,d) PRODUCTION over the full proton-target run set -> reco_gf/.
# 4-parallel (CPU-bound fits, ~0.8 GB each). Resumable (skips _genfit.root >1MB).
# Args: $1=minIter (default 2)  $2=maxIter (default 5)  -- set from the speed scan.
MI=${1:-2}; MA=${2:-5}
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF
source build/config.sh >/dev/null 2>&1
export GENFIT=/home/yassid/fair_install/GenFit
export LD_LIBRARY_PATH=$GENFIT/lib:$LD_LIBRARY_PATH
cd macro/Unpack_HDF5/a1975/UKF
mkdir -p /mnt/f/a1975/reco_gf
export RECO=/mnt/f/a1975/reco/ OUT=/mnt/f/a1975/reco_gf/ MI MA
fit_one() {
  r=$1; run=run_$r; outf=${OUT}${run}_genfit.root
  [ -f "${RECO}${run}_reco.root" ] || { echo ">>> $run no reco, skip"; return; }
  if [ -f "$outf" ] && [ $(stat -c %s "$outf") -gt 1000000 ]; then echo ">>> $run done, skip"; return; fi
  echo ">>> $run (genfit iter $MI-$MA) ..."
  root -b -q "pipeline/fitGenfit_a1975.C(\"$run\", -1, \"deuteron\", 2.85, 0.083147, \"$RECO\", \"\", \"$OUT\", \"deuteron_H2_catima.txt\", $MI, $MA)" \
    > /tmp/genfitprod_${r}.log 2>&1
  echo "    done $run ($(grep -iE 'Real' /tmp/genfitprod_${r}.log | tail -1))"
}
export -f fit_one
# full proton-target set: runs 0106-0189
RUNS=$(seq -w 106 189)
printf '%s\n' $RUNS | sed 's/^/0/' | xargs -P4 -I{} bash -c 'fit_one {}'
echo "=== GENFIT PRODUCTION DONE (iter $MI-$MA) ==="
