#!/usr/bin/env bash
# Open the a1975 (p,p) proton gate editor.
#
#   ./pid/draw_gate.sh                      # runs 106-130, writes pid/proton_hand.json
#   ./pid/draw_gate.sh 106 189 pid/mine.json
#
# Run it FROM macro/Unpack_HDF5/a1975/UKF -- the gate paths inside the macro, and the ones
# fitGenfitter_a1975.C resolves later, are relative to that directory.
#
# In the window: [Draw new gate] -> click the polygon vertices -> DOUBLE-CLICK to close ->
# [Evaluate] -> [Save JSON]. Evaluate reports, from the simulation, what fraction of REAL protons
# the polygon keeps, how many data tracks it admits against the current production gate, and how
# much of that overlaps the deuteron gate. Draw, evaluate, redraw until those three read well.
#
# Once saved, apply it with:  ./pp/refit_pp.sh pid/proton_hand.json _pphand 4
set -uo pipefail
LO=${1:-106}; HI=${2:-130}; OUT=${3:-pid/proton_hand.json}
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
UKF=$REPO/macro/Unpack_HDF5/a1975/UKF
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
cd "$UKF"
[ -n "${DISPLAY:-}" ] || echo "WARNING: DISPLAY is unset; the GUI needs WSLg. ROOT TGMainFrame works there, Eve does not."
exec root -l "pid/gate_draw_pp.C($LO,$HI,\"$OUT\")"
