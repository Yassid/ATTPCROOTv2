#!/usr/bin/env bash
# UKF proton refit of the 47 D2 (d,p) recos (reuses existing multifit recos -> fast), 4 cores,
# resumable via .ukfdone. Then ex_dp with the UKF fits, saved separately for genfit-vs-UKF compare.
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
OUT=/mnt/f/a1975/reco_d2/
LOGD=/tmp/claude-1000/-home-yassid/b5b91330-d765-4197-9f1c-6da8c3d5a215/scratchpad/ukf
mkdir -p "$LOGD"
cd "$REPO/macro/Unpack_HDF5/a1975/D2_UKF"
source "$REPO/build/config.sh" >/dev/null 2>&1
NUMS="0016 0017 0018 0019 0020 0021 0022 0023 0026 0027 0031 0032 0034 0036 0037 0038 0039 0040 0041 0042 0043 0044 0046 0048 0057 0058 0076 0077 0078 0079 0080 0082 0083 0084 0085 0086 0087 0088 0089 0091 0092 0095 0096 0097 0098 0102 0103"
NPAR=4
ukffit(){
  local r=run_$1_multifit
  local o=${OUT}${r}_genfitter_p_UKF.root
  [ -f "${o}.ukfdone" ] && return
  [ -f "${OUT}${r}_reco.root.mfdone" ] || { echo "$(date +%H:%M) NO-RECO $1"; return; }
  root -l -b -q "fitUKF_a1975_deuterium.C(\"$r\",-1,\"proton\",-1,2.85,9e-5,\"$OUT\")" > "$LOGD/$r.log" 2>&1
  grep -qi "Done\." "$LOGD/$r.log" && touch "${o}.ukfdone" || echo "$(date +%H:%M) UKF-FAIL $1"
}
echo "=== UKF proton refit: 47 runs, ${NPAR}-par $(date) ==="
cp -f proton_kin_dp.root proton_kin_dp_genfit.root 2>/dev/null   # preserve the genfit cache
for n in $NUMS; do ukffit "$n" & while [ "$(jobs -r | wc -l)" -ge $NPAR ]; do sleep 3; done; done
wait
ndone=$(ls ${OUT}*_multifit_genfitter_p_UKF.root.ukfdone 2>/dev/null | wc -l)
echo "=== UKF fits done: ${ndone}/47. ex_dp (UKF) $(date) ==="
RUNS_MF=""; for n in $NUMS; do [ -f "${OUT}run_${n}_multifit_genfitter_p_UKF.root.ukfdone" ] && RUNS_MF="$RUNS_MF,run_${n}_multifit"; done
RUNS_MF=${RUNS_MF#,}
root -l -b -q "ex_dp_a1975.C(\"$RUNS_MF\",\"$OUT\",\"_UKF\")" > "$LOGD/ex_dp_ukf.log" 2>&1
grep -iE "candidates|17C Ex:" "$LOGD/ex_dp_ukf.log" | tail -1
cp -f proton_kin_dp.root proton_kin_dp_UKF.root
cp -f plots/ex_dp_spectrum.png plots/ex_dp_spectrum_UKF.png
cp -f proton_kin_dp_genfit.root proton_kin_dp.root       # restore genfit as the canonical cache
echo "UKF_REFIT_DONE $(date)  -> proton_kin_dp_UKF.root + plots/ex_dp_spectrum_UKF.png"
