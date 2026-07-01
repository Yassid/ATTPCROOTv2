#!/usr/bin/env bash
# Overnight FULL D2 (d,p)17C production with the CURRENT pipeline:
#   AtPSAMultiFit (thr40, full cloud) -> HDBSCAN mover (+ clusterize fix) -> genfit proton B=-2.85
# 47 canonical (d,p) runs, ALL events. Per-run reco+fit job, NPAR concurrent, RESUMABLE via
# .mfdone markers. Logs -> scratchpad (keeps the folder clean). Then ex_dp -> 17C Ex over all runs.
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
OUT=/mnt/f/a1975/reco_d2/
H5=/mnt/f/a1975/h5/
LOGD=/tmp/claude-1000/-home-yassid/b5b91330-d765-4197-9f1c-6da8c3d5a215/scratchpad/dp_full
mkdir -p "$LOGD"
cd "$REPO/macro/Unpack_HDF5/a1975/D2_UKF"
source "$REPO/build/config.sh" >/dev/null 2>&1
NUMS="0016 0017 0018 0019 0020 0021 0022 0023 0026 0027 0031 0032 0034 0036 0037 0038 0039 0040 0041 0042 0043 0044 0046 0048 0057 0058 0076 0077 0078 0079 0080 0082 0083 0084 0085 0086 0087 0088 0089 0091 0092 0095 0096 0097 0098 0102 0103"
NPAR=3

recofit(){
  local r=run_$1
  local reco=${OUT}${r}_multifit_reco.root
  local fit=${OUT}${r}_multifit_genfitter_p.root
  if [ ! -f "${reco}.mfdone" ]; then
    root -l -b -q "unpackReco_multifit.C(\"$r\",0,false,\"$OUT\",\"$H5\",false,true,\"multifit\",0,40,\"hdbscan\",0,3)" > "$LOGD/reco_$r.log" 2>&1
    grep -q "Done ->" "$LOGD/reco_$r.log" && touch "${reco}.mfdone" || { echo "$(date +%H:%M) RECO-FAIL $r"; return; }
  fi
  if [ ! -f "${fit}.mfdone" ]; then
    root -l -b -q "fitGenfitter_a1975_deuterium.C(\"${r}_multifit\",-1,\"$OUT\",\"\",\"$OUT\",-2.85)" > "$LOGD/fit_$r.log" 2>&1
    grep -qi "Done\." "$LOGD/fit_$r.log" && touch "${fit}.mfdone" || echo "$(date +%H:%M) FIT-FAIL $r"
  fi
  echo "$(date +%H:%M) DONE $r"
}

echo "=== OVERNIGHT FULL D2 (d,p): 47 runs, ${NPAR}-par, multifit/HDBSCAN/B=-2.85 $(date) ==="
for n in $NUMS; do
  recofit "$n" &
  while [ "$(jobs -r | wc -l)" -ge $NPAR ]; do sleep 5; done
done
wait
nfit=$(ls ${OUT}*_multifit_genfitter_p.root.mfdone 2>/dev/null | wc -l)
echo "=== reco+fit complete: ${nfit}/47 fits done. ex_dp $(date) ==="
RUNS_MF=""; for n in $NUMS; do [ -f "${OUT}run_${n}_multifit_genfitter_p.root.mfdone" ] && RUNS_MF="$RUNS_MF,run_${n}_multifit"; done
RUNS_MF=${RUNS_MF#,}
root -l -b -q "ex_dp_a1975.C(\"$RUNS_MF\",\"$OUT\")" > "$LOGD/ex_dp.log" 2>&1
grep -iE "candidates|17C Ex:|saved" "$LOGD/ex_dp.log" | tail -3
echo "OVERNIGHT_FULL_DONE $(date)"
