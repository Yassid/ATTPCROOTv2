#!/usr/bin/env bash
# SENSITIVITY CHECK: does the assumed per-hit measurement error hide the gain from 2 mm pads?
#
#   ./meassigma_check.sh [state] [B_T] [pad_mm] [seed] [nEvents] [sigma_list]
#   ./meassigma_check.sh gs 7.0 2.0 9301 4000 "4.0 2.0 1.0"
#
# The campaign holds genfit's measSigma at 4.0 mm in every cell, because it is an analysis
# parameter and the campaign varies the field and the pitch only. But 4 mm was chosen for
# 8 x 12 mm pads with 3.7 mm of transverse diffusion; at 2 mm pitch and 7 T the true hit error
# is closer to 1 mm, and a Kalman filter told its measurements are four times worse than they
# are leans on the model instead of the data. If the fine-pitch resolution improves when
# measSigma is lowered, then the campaign numbers for those cells are a LOWER BOUND on the gain
# and should be quoted that way.
#
# Runs its own generation and reconstruction once, then only re-fits -- the reco is the
# expensive stage and it does not depend on measSigma.
set -eo pipefail
STATE=${1:-gs}; BT=${2:-7.0}; PAD=${3:-2.0}; SEED=${4:-9301}; NEV=${5:-4000}
SIGMAS=${6:-"4.0 2.0 1.0"}
case "$STATE" in
   gs) EX=0.0 ;; ex6094) EX=6.094 ;; ex6728) EX=6.728 ;; ex7012) EX=7.012 ;; ex8317) EX=8.317 ;;
   *) echo "unknown state '$STATE'"; exit 2 ;;
esac
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
SIM="$REPO/macro/Simulation/ATTPC/14C_pp"; HF="$SIM/highfield"
UKF="$REPO/macro/Unpack_HDF5/a1954/UKF"
BTAG=$(awk -v b="$BT" 'BEGIN{printf "b%03d", b*100}')
OUT=${MS_OUT:-/mnt/f/a1954_C14_hf_meassigma}; mkdir -p "$OUT/work/data"
BKG=$(awk -v b="$BT" 'BEGIN{printf "%.1f", -10*b}')
BNEG=$(awk -v b="$BT" 'BEGIN{printf "%.2f", -b}')
PAR="ATTPC.a1954_C14_hf_${BTAG}.par"
J="${STATE}_s${SEED}_${BTAG}"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"

[ -s "$OUT/${J}_sim.root" ] || ( cd "$OUT/work" && root -b -q -l "$SIM/C14_pp_sim.C($NEV,2.0,178.0,\"TGeant4\",$BKG,\"$OUT/${J}_sim.root\",$EX,$SEED,\"ATTPC_H300torr_RT.root\")" ) > "$OUT/${J}_gen.log" 2>&1
[ -s "$OUT/${J}_reco.root" ] || ( cd "$SIM" && root -b -q -l "run_reco_C14.C(\"$OUT/${J}_sim.root\",\"$OUT/${J}_reco.root\",\"$PAR\",20,20,8,0,30.0,\"mover\",$PAD,500.0,kFALSE)" ) > "$OUT/${J}_reco.log" 2>&1
[ -s "$OUT/${J}_reco.root" ] || { echo "RECO_FAILED"; exit 1; }

cd "$SIM"
for MS in $SIGMAS; do
   T="${J}_ms${MS}"
   # fitGenfit_C14.C names its output <fileName>_genfit.root, so the input link has to carry the
   # per-sigma tag as well or the three fits overwrite each other
   ln -sf "$OUT/${J}_reco.root" "$OUT/${T}_reco.root"
   [ -s "$OUT/${T}_genfit.root" ] || root -b -q -l "$UKF/pipeline/fitGenfit_C14.C(\"$T\",-1,\"$OUT/\",\"\",\"$OUT/\",$BNEG,2,5,\"\",$MS,10.0,170.0,kTRUE,kFALSE,\"proton\",\"ATTPC_H300torr_RT\",kFALSE,kTRUE,0.0,1,kTRUE,kFALSE,kFALSE,kFALSE)" > "$OUT/${T}_fit.log" 2>&1
   # chi2/ndf scales as 1/measSigma^2, so the cut cannot be held at 5 across the scan without
   # meaning something different in each arm. Quote both: the fixed cut and no cut at all.
   echo "===== measSigma = $MS mm, chi2/ndf < 5 ====="
   root -b -q -l "$HF/ex_res_C14_hf.C(\"$OUT/${J}_sim.root\",\"$OUT/${T}_genfit.root\",\"${T}_c5\",$EX,159.75,5.0,kTRUE,\"$OUT/\")" 2>&1 | sed -n '/Ex resolution/,/^$/p'
   echo "===== measSigma = $MS mm, no chi2 cut ====="
   root -b -q -l "$HF/ex_res_C14_hf.C(\"$OUT/${J}_sim.root\",\"$OUT/${T}_genfit.root\",\"${T}_nc\",$EX,159.75,1e9,kTRUE,\"$OUT/\")" 2>&1 | sed -n '/Ex resolution/,/^$/p'
done
echo "meassigma check done -> $OUT"
