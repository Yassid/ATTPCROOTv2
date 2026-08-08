#!/usr/bin/env bash
# Second half of the triplclust-vs-HDBSCAN test: wait for pra_cmp_C14.sh to finish both recos,
# fit each with the SAME UKF configuration, build an Ex cache for each, and print the elastic
# ridge against the kinematic line so the two clusterers can be compared directly.
#
# The comparison that matters is "ridge - elastic line" as a function of theta_lab: that is the
# +7-19 % bias for theta_lab < 62 deg that neither the beam energy, the drift velocity, the field
# scale nor the simulation reproduces. If triplclust removes it, the clusterer is the cause.
#
#   ./pra_cmp_finish_C14.sh run_0058

set -uo pipefail
RUN="${1:-run_0058}"
BASE="/home/yassid/a1954_C14_pracmp"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"

set +u
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
set -u

# wait for BOTH markers (the marker means completed, not merely "file exists")
echo "waiting for reco markers..."
for i in $(seq 1 720); do
   if [[ -f "$BASE/tc/${RUN}.marker" && -f "$BASE/hdbscan/${RUN}.marker" ]]; then
      echo "both recos done"
      break
   fi
   sleep 30
done
if [[ ! -f "$BASE/tc/${RUN}.marker" || ! -f "$BASE/hdbscan/${RUN}.marker" ]]; then
   echo "FATAL: timed out waiting for reco markers" >&2
   exit 1
fi

for pra in tc hdbscan; do
   d="$BASE/$pra"
   echo "=== fitting $RUN ($pra) ==="
   root -b -q -l "$HERE/pipeline/fitUKF_C14.C(\"$RUN\",-1,\"proton\",-1,2.85,3.553e-5,\"\",\"$d/\",0.5,0.1,1,10,\"$d/\")" \
      > "$d/logs/${RUN}_fit.log" 2>&1
   if [[ ! -s "$d/${RUN}_ukf.root" ]]; then
      echo "FATAL: no ${RUN}_ukf.root for $pra -- see $d/logs/${RUN}_fit.log" >&2
      exit 1
   fi
   echo "=== Ex cache ($pra) ==="
   root -b -q -l "$HERE/pp/ex_C14.C(\"$RUN\",\"$d/\",161.0,1e9,\"_pracmp_${pra}\",1.007825,14.003242,\"\",\"ukf\",kFALSE)" \
      > "$d/logs/${RUN}_ex.log" 2>&1
   grep -E "good track|protons" "$d/logs/${RUN}_ex.log" | tail -1
done

echo "=== ridge vs kinematic line, tc vs hdbscan ==="
root -b -q -l "$HERE/pp/pra_ridge_C14.C()" 2>&1 | sed -n '/=====/,$p'
echo "PRA COMPARISON DONE"
