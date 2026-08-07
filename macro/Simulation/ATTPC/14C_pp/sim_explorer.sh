#!/bin/bash
# Build the browser Ex explorer for the 14C(p,p') SIMULATION, the sim twin of the data page
# produced by macro/Unpack_HDF5/a1954/UKF/pp/make_explorer_html.C.
#
#   ./sim_explorer.sh [NAME] [OUT]
#     NAME  fit basename in ./data/  (default simg = the TRUTH-GATED protons; sim = ungated)
#     OUT   output html               (default ~/a1954_C14_sim_explorer.html)
#
# Differences from the data page, both deliberate:
#   * chi2Cut = 1e9 when the cache is built, so the page's own chi2/ndf slider spans the whole
#     sample. The data caches were built at chi2 < 5, which freezes that control.
#   * reference levels = truth only. C14_pp_sim.C generates ELASTIC scattering (ExE = 0 for
#     every state), so overlaying the 6.09/6.59/7.01 MeV 14C levels of the data page would
#     draw loci that nothing in this sample can populate.
#
# The closure number to read off the page is the g.s. centroid: truth is EXACTLY Ex = 0 at
# Ebeam = 161.00 MeV, so mu IS the chain bias.
set -eo pipefail

NAME=${1:-simg}
OUT=${2:-$HOME/a1954_C14_sim_explorer.html}

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../../../.." && pwd)
UKFDIR="$REPO/macro/Unpack_HDF5/a1954/UKF"
DATA="$HERE/data/"

set +u
source "$REPO/build/config.sh" >/dev/null 2>&1   # set -u + thisroot.sh kills the script
set -u

for F in ukf genfit; do
  [ -f "$DATA/${NAME}_$F.root" ] || { echo "MISSING $DATA/${NAME}_$F.root -- run fit_sim_C14.sh first"; exit 1; }
done

cd "$UKFDIR"
for F in ukf genfit; do
  echo "[$(date +%H:%M:%S)] cache: ${NAME}_$F  (no chi2 cut)"
  root -b -q -l "pp/ex_C14.C(\"$NAME\",\"$DATA\",161.0,1e9,\"_${NAME}x_$F\",1.007825,14.003242,\"14C(p,p') simulation\",\"$F\")" 2>&1 \
    | grep -E "good track|saved"
done

echo "[$(date +%H:%M:%S)] baking $OUT"
root -b -q -l "pp/make_explorer_html.C(\"pp/plots/proton_kin_${NAME}x_ukf.root\",\"$OUT\",\"14C(p,p') SIMULATION\",161.0,14.003242,1.007825,1.007825,14.003242,14,\"0:truth g.s.\",\"pp/plots/proton_kin_${NAME}x_genfit.root\")" 2>&1 \
  | grep -E "wrote|UKF|GENFIT"
