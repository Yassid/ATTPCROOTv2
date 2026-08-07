#!/usr/bin/env bash
# Full chain with the CORRECTED field convention (B = -28.5 kG, i.e. the data's sign):
#   generate -> reco (3 cm hole, mover join) -> truth gate -> clustering score -> PID
# Marker written only after every stage has produced real output.
set -eo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); SIM=$(cd "$HERE/.." && pwd)
REPO=$(cd "$SIM/../../../.." && pwd); OUT="$HERE/negB"; mkdir -p "$OUT"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
cd "$SIM"
echo "[$(date +%H:%M:%S)] generating 8000 events at B = -28.5 kG"
root -b -q -l "C14_pp_sim.C(8000,5.0,120.0,\"TGeant4\",-28.5,\"$OUT/attpcsim_negB.root\")" > "$OUT/gen.log" 2>&1
[ -s "$OUT/attpcsim_negB.root" ] || { echo GEN_EMPTY; exit 1; }
echo "[$(date +%H:%M:%S)] reco (hole 30 mm, join mover)"
root -b -q -l "run_reco_C14.C(\"$OUT/attpcsim_negB.root\",\"$OUT/negB_reco.root\",\"ATTPC.a1954_C14_sim.par\",20,20,8,0,30.0,\"mover\")" > "$OUT/reco.log" 2>&1
[ -s "$OUT/negB_reco.root" ] || { echo RECO_EMPTY; exit 1; }
echo "[$(date +%H:%M:%S)] truth gate + clustering score"
root -b -q -l "gate_truth_C14.C(\"$OUT/negB_reco.root\",\"$OUT/\",\"negBg\",0.6,4)" > "$OUT/gate.log" 2>&1
root -b -q -l "cluster_eval_C14.C(\"$OUT/negB_reco.root\",\"negB\",-1)" > "$OUT/clus.log" 2>&1
grep -qE "RECOIL PROTON" "$OUT/clus.log" || { echo SCORE_FAILED; exit 1; }
echo "[$(date +%H:%M:%S)] PID"
cd "$REPO/macro/Unpack_HDF5/a1954/UKF"
root -b -q -l "pid/pid_C14.C(\"negBg\",\"$OUT/\",-1,2.85,\"_negBg\",\"pid/proton_14C.json\",80,2.5,4000)" > "$OUT/pid.log" 2>&1
grep -qE "tracks=" "$OUT/pid.log" || { echo PID_FAILED; exit 1; }
echo COMPLETED > "$OUT/negB.marker"
echo "[$(date +%H:%M:%S)] done"
