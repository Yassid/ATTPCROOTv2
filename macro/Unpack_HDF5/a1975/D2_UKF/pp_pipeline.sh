#!/usr/bin/env bash
# 16C+p (proton-target) with the NEW setup: multifit thr40 + HDBSCAN mover reco (ATTPC.a1954.par)
# -> genfit proton -> (p,p) elastic kinematics. Small run set to demonstrate; resumable via .ppdone.
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
OUT=/mnt/f/a1975/reco_pp/
H5=/mnt/f/a1975/h5/
PAR="ATTPC.a1954.par"
LOGD=/tmp/claude-1000/-home-yassid/b5b91330-d765-4197-9f1c-6da8c3d5a215/scratchpad/pp
mkdir -p "$OUT" "$LOGD"
cd "$REPO/macro/Unpack_HDF5/a1975/D2_UKF"
source "$REPO/build/config.sh" >/dev/null 2>&1
NUMS="0106 0107 0108 0109 0110"
NPAR=2
recofit(){
  local r=run_$1
  local reco=${OUT}${r}_multifit_reco.root
  local fit=${OUT}${r}_multifit_genfitter_p.root
  if [ ! -f "${reco}.ppdone" ]; then
    root -l -b -q "unpackReco_multifit.C(\"$r\",0,false,\"$OUT\",\"$H5\",false,true,\"multifit\",0,40,\"hdbscan\",0,3,0,0.1,\"$PAR\")" > "$LOGD/reco_$r.log" 2>&1
    grep -q "Done ->" "$LOGD/reco_$r.log" && touch "${reco}.ppdone" || { echo "$(date +%H:%M) RECO-FAIL $r"; return; }
  fi
  if [ ! -f "${fit}.ppdone" ]; then
    root -l -b -q "fitGenfitter_a1975_deuterium.C(\"${r}_multifit\",-1,\"$OUT\",\"\",\"$OUT\",-2.85)" > "$LOGD/fit_$r.log" 2>&1
    grep -qi "Done\." "$LOGD/fit_$r.log" && touch "${fit}.ppdone" || echo "$(date +%H:%M) FIT-FAIL $r"
  fi
  echo "$(date +%H:%M) DONE $r"
}
echo "=== 16C+p (p,p) NEW setup: 5 runs, ${NPAR}-par, a1954.par $(date) ==="
for n in $NUMS; do recofit "$n" & while [ "$(jobs -r | wc -l)" -ge $NPAR ]; do sleep 5; done; done
wait
RUNS=""; for n in $NUMS; do [ -f "${OUT}run_${n}_multifit_genfitter_p.root.ppdone" ] && RUNS="$RUNS,run_${n}_multifit"; done
RUNS=${RUNS#,}
echo "=== ex_elastic_pp on: $RUNS $(date) ==="
root -l -b -q "ex_elastic_pp.C(\"$RUNS\",\"$OUT\")" > "$LOGD/ex_pp.log" 2>&1
grep -iE "candidates|Ex:|saved" "$LOGD/ex_pp.log" | tail -3
echo "PP_PIPELINE_DONE $(date)"
