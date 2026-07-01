#!/usr/bin/env bash
# Deuteron-hypothesis genfit of the D2 recos (for 16C(d,d) ELASTIC), 2 cores (leaves 4 for UKF),
# resumable via .ddone. Then ex_elastic_dd with the deuteron gate -> 16C Ex (=0 for elastic).
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
OUT=/mnt/f/a1975/reco_d2/
LOGD=/tmp/claude-1000/-home-yassid/b5b91330-d765-4197-9f1c-6da8c3d5a215/scratchpad/deut
mkdir -p "$LOGD"
cd "$REPO/macro/Unpack_HDF5/a1975/D2_UKF"
source "$REPO/build/config.sh" >/dev/null 2>&1
NUMS="0016 0017 0018 0019 0020 0021 0022 0023 0026 0027 0031 0032 0034 0036 0037 0038 0039 0040 0041 0042 0043 0044 0046 0048 0057 0058 0076 0077 0078 0079 0080 0082 0083 0084 0085 0086 0087 0088 0089 0091 0092 0095 0096 0097 0098 0102 0103"
NPAR=2
dfit(){
  local r=run_$1_multifit
  local o=${OUT}${r}_genfitter_d.root
  [ -f "${o}.ddone" ] && return
  [ -f "${OUT}${r}_reco.root.mfdone" ] || return
  root -l -b -q "fitGenfitter_a1975_deuteron.C(\"$r\",-1,\"$OUT\",\"\",\"$OUT\",-2.85)" > "$LOGD/$r.log" 2>&1
  grep -qi "Done\." "$LOGD/$r.log" && touch "${o}.ddone" || echo "$(date +%H:%M) DFIT-FAIL $1"
}
echo "=== deuteron fit (elastic): 47 runs, ${NPAR}-par $(date) ==="
for n in $NUMS; do dfit "$n" & while [ "$(jobs -r | wc -l)" -ge $NPAR ]; do sleep 4; done; done
wait
nd=$(ls ${OUT}*_multifit_genfitter_d.root.ddone 2>/dev/null | wc -l)
echo "=== deuteron fits done: ${nd}/47. ex_elastic_dd $(date) ==="
RUNS=""; for n in $NUMS; do [ -f "${OUT}run_${n}_multifit_genfitter_d.root.ddone" ] && RUNS="$RUNS,run_${n}_multifit"; done
RUNS=${RUNS#,}
root -l -b -q "ex_elastic_dd.C(\"$RUNS\",\"$OUT\")" > "$LOGD/ex_elastic.log" 2>&1
grep -iE "candidates|Ex:|saved" "$LOGD/ex_elastic.log" | tail -3
echo "DEUTERON_ELASTIC_DONE $(date)"
