#!/bin/bash
# 14C(p,t)12C browser explorer, on the 14-run GENFIT+CATIMA triton production.
#
# THE LEVEL SCHEME MUST BE PASSED EXPLICITLY. The generator infers it from the residual mass and
# its fallback is 13C -- that is how 12Be once got 13C lines drawn on it. The residual here is 12C.
#
# KINEMATIC CEILING: Q(p,t) = -4.641 MeV and Ecm = 10.72 MeV, so only 6.08 MeV of excitation is
# available. The 7.654 and 9.641 lines are drawn anyway, deliberately: they are CLOSED, so any
# yield appearing on them is a defect and the lines are there to make that visible.
REPO="/home/yassid/fair_install/ATTPCROOTv2-OpenKF"
PL="$REPO/macro/Unpack_HDF5/a1954/UKF/pp/plots"
V="/home/yassid/a1975_analysis/common/macros/make_explorer_html.C"
OUT="${1:-/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/pages/a1954_C14_pt_explorer.html}"
mkdir -p "$(dirname "$OUT")"
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
# 12C levels. 0 and 4.4398 are open; 7.6542 and 9.641 are above the ceiling and must stay empty.
REFEX="0,4.4398,7.6542,9.641"
root -b -l -q "$V(\"$PL/proton_kinpt14.root\",\"$OUT\",\"14C(p,t)12C\",159.75,14.003242,1.007825,3.016049,12.000000,14,\"$REFEX\",\"\",\"\",\"genfit+CATIMA, 14 runs\")" 2>&1 | grep -E "tracks|open it|wrote"
python3 /home/yassid/a1975_analysis/common/macros/add_keoff.py "$OUT" 2>&1 | grep -E "patched|WARNING"
# the (p,t) Ex range is -3..11, not the (p,p') -5..25
sed -i 's|id="exLo" step="0.5" value="-5"|id="exLo" step="0.5" value="-3"|' "$OUT" 2>/dev/null
sed -i 's|id="exHi" step="0.5" value="25"|id="exHi" step="0.5" value="11"|' "$OUT" 2>/dev/null
cp "$OUT" "/mnt/c/Users/Yassid/Desktop/" 2>/dev/null && echo "-> Desktop: $(basename $OUT)"
echo "$OUT"
