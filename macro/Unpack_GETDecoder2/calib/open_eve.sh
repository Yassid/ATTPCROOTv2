#!/bin/bash
# Open the TEve display on a Dec2014 alpha run (OpenKF-Claude branch).
# Needs a tty -- launch inside a terminal window, not backgrounded.
#   usage: open_eve.sh [file.root]
F=${1:-alpha_run_0080.root}
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh
source /home/yassid/fair_install/ATTPCROOTv2/build/config.sh
cd /home/yassid/fair_install/ATTPCROOTv2/macro/Unpack_GRAW_Dec2014 || exit 1
echo "=== opening $F in the Eve viewer ==="
root -l "eve_Dec2014_alphas.C(\"$F\")"
