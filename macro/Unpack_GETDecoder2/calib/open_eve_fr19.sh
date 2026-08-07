#!/bin/bash
# Open the PCL-free TEve display on a Dec2014 alpha file unpacked with the fr19port chain.
# Needs a tty -- launch inside a terminal window, not backgrounded.
#   usage: open_eve_fr19.sh [file.root]
F=${1:-/home/yassid/dec2014_alphas_reco/lowP/alpha_run_0128_hits.root}
source /home/yassid/fair_install/ATTPCROOTv2_fr19port/setup_fr19port.sh > /dev/null 2>&1
cd /home/yassid/fair_install/ATTPCROOTv2_fr19port/macro/Unpack_GETDecoder2 || exit 1
echo "=== opening $F in the fr19port Eve viewer ==="
root -l "run_eve_Dec2014_alphas.C(\"$F\")"
