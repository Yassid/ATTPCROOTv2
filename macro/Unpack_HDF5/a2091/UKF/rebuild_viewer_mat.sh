#!/bin/bash
# Wait for the genfit matEffects=kTRUE pass to finish, then rebuild the 300-torr explorer
# with the material-corrected GENFIT cache and stage it for the Windows browser.
# Absolute paths everywhere -- a relative ./ launch is what silently killed the previous
# chained job when the caller's cwd had moved.
set -u
REPO="/home/yassid/fair_install/ATTPCROOTv2"
HERE="$REPO/macro/Unpack_HDF5/a2091/UKF"
P="$HERE/pp/plots"
MATDIR="/home/yassid/a2091_C15_genfit_mat/"
CHAINLOG="/home/yassid/a2091_C15_genfit_mat_chain.log"
OUT="/home/yassid/a2091_C15_pp_300torr.html"
WINHOME="/mnt/c/Users/Yassid"
RUNS="run_0138,run_0056,run_0057,run_0058,run_0059,run_0060,run_0061,run_0062,run_0063,run_0064,run_0065,run_0066,run_0068,run_0069"

echo "[$(date +%H:%M:%S)] waiting for the genfit matEffects pass..."
while ! grep -q "GENFIT matEffects PASS DONE" "$CHAINLOG" 2>/dev/null; do sleep 15; done
N=$(ls "$MATDIR"*_genfit.root 2>/dev/null | wc -l)
echo "[$(date +%H:%M:%S)] pass done, $N genfit files present"

source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
cd "$HERE" || exit 1

echo "[$(date +%H:%M:%S)] building material-corrected GENFIT Ex cache..."
root -b -q -l "pp/ex_C15.C(\"$RUNS\",\"$MATDIR\",195.0,5.0,\"_300_genfit_mat\",1.007825,15.0105993,\"\",\"genfit\")" 2>&1 | grep "good track"

echo "[$(date +%H:%M:%S)] rebuilding the explorer..."
root -b -q -l "pp/make_explorer_html.C(\"$P/proton_kin_300_ukf.root\",\"$OUT\",\"15C(p,p') 300 torr + matEffects - 14 runs\",161,15.0105993,1.007825,1.007825,15.0105993,14,\"6.09,6.59,6.73,6.90,7.01,7.34\",\"$P/proton_kin_300_genfit_mat.root\")" 2>&1 | grep -E "wrote|UKF |GENFIT "

cp "$OUT" "$WINHOME/" 2>/dev/null
cp "$OUT" "$WINHOME/Desktop/" 2>/dev/null

echo "[$(date +%H:%M:%S)] ===== headline numbers ====="
for c in 300_ukf 300_genfit 300_genfit_mat; do
  printf "  %-16s " "$c"
  root -b -q -l "pp/hyp_C15.C(\"plots/proton_kin_$c.root\",161,\"$c\")" 2>&1 | grep "(0) as measured" | sed 's/(0) as measured *: //'
done
echo "[$(date +%H:%M:%S)] VIEWER REBUILD DONE -> $OUT"
