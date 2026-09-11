#!/bin/bash
# Dec2014 alpha kinematics drawer -- ONE window.
#   ./open_kin.sh [vD] [pressure] [Ebeam]
#
# The launch recipe, and why each half is needed:
#   * the MACRO calls gApplication->Run() itself.  Passing -e 'gApplication->Run()' instead
#     fails with "__clang_call_terminate unresolved".
#   * stdin must be HELD OPEN (tail -f /dev/null).  With stdin at /dev/null, ROOT sees EOF,
#     Run() returns immediately, and the macro is re-processed -- a new window every ~0.4 s.
#   * /tmp/kin_draw.lock is a belt-and-braces guard inside the macro: if anything ever
#     respawns again, the second instance refuses to start instead of flooding the screen.
cd "$(dirname "$0")"
rm -f /tmp/kin_draw.lock
trap 'rm -f /tmp/kin_draw.lock' EXIT
VD=${1:-2.142}; P=${2:-299.5}; EB=${3:-11.40}
source /home/yassid/fair_install/ATTPCROOTv2_fr19port/setup_fr19port.sh >/dev/null 2>&1
[ -s kin_points.root ] || { echo "kin_points.root missing"; exit 1; }
echo "launching: v_D=$VD, $P torr, E_beam=$EB   (Quit button to close)"
tail -f /dev/null | root -l "kin_draw.C(\"kin_points.root\",\"/home/yassid/dec2014_kin/eloss_alpha_heco2_11p4.txt\",$VD,$P,$EB)"
