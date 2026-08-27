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
# TWO SETS: the adopted [5,175] fit and the [10,170] baseline, so the theta-window change can be
# flipped between in the browser instead of being taken on trust. The proton window [10,170] slices
# the triton band (reconstructed 156-172) and costs 12 % of it, all at the forward end.
root -b -l -q "$V(\"$PL/proton_kinpt14t.root\",\"$OUT\",\"14C(p,t)12C\",159.75,14.003242,1.007825,3.016049,12.000000,14,\"$REFEX\",\"$PL/proton_kinpt14w.root\",\"$PL/proton_kinpt14.root\",\"tight gate th[5|175],broad gate th[5|175],broad gate th[10|170]\")"
python3 /home/yassid/a1975_analysis/common/macros/add_keoff.py "$OUT" 2>&1 | grep -E "patched|WARNING"
# The generator hardcodes its own experiment in the header eyebrow -- it ships "a2091", which is a
# different experiment entirely. Left alone the page announces itself as a2091 data.
sed -i 's/"eyebrow":"a2091 . AT-TPC"/"eyebrow":"a1954 . AT-TPC"/' "$OUT"
# the (p,t) Ex range is -3..11, not the (p,p') -5..25
sed -i 's|id="exLo" step="0.5" value="-5"|id="exLo" step="0.5" value="-6"|' "$OUT" 2>/dev/null
sed -i 's|id="exHi" step="0.5" value="25"|id="exHi" step="0.5" value="16"|' "$OUT" 2>/dev/null
cp "$OUT" "/mnt/c/Users/Yassid/Desktop/" 2>/dev/null && echo "-> Desktop: $(basename $OUT)"
echo "$OUT"
