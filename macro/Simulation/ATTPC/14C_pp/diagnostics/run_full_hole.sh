#!/usr/bin/env bash
# Full 8000-event sim reco with the 3 cm beam hole, WITH and WITHOUT the HDBSCAN join stage.
# Markers are written ONLY after the whole chain (reco -> truth gate -> cluster score) succeeds,
# so "marker exists" means COMPLETED, never "a file appeared".
set -eo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); SIM=$(cd "$HERE/.." && pwd)
REPO=$(cd "$SIM/../../../.." && pwd); OUT="$HERE/full"; mkdir -p "$OUT"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
TAG=$1; JOIN=$2
cd "$SIM"
root -b -q -l "run_reco_C14.C(\"./data/attpcsim.root\",\"$OUT/${TAG}_reco.root\",\"ATTPC.a1954_C14_sim.par\",20,20,8,0,30.0,\"$JOIN\")" > "$OUT/${TAG}_reco.log" 2>&1
[ -s "$OUT/${TAG}_reco.root" ] || { echo "$TAG RECO_EMPTY"; exit 1; }
root -b -q -l "gate_truth_C14.C(\"$OUT/${TAG}_reco.root\",\"$OUT/\",\"${TAG}g\",0.6,4)" > "$OUT/${TAG}_gate.log" 2>&1
root -b -q -l "cluster_eval_C14.C(\"$OUT/${TAG}_reco.root\",\"$TAG\",-1)" > "$OUT/${TAG}_clus.log" 2>&1
grep -qE "RECOIL PROTON" "$OUT/${TAG}_clus.log" || { echo "$TAG SCORE_FAILED"; exit 1; }
echo COMPLETED > "$OUT/${TAG}.marker"
