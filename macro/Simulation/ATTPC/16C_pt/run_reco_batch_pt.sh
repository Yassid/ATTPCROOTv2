#!/usr/bin/env bash
# Reconstruct the five 16C(p,t)14C simulated states at the ADOPTED gain.
#
# GAIN 25000, MEASURED not inherited (make_gain_pars_pt.sh + hits_per_track.C, against H2 runs
# 0106/0107/0120 whose median hits per track is 45/48/53):
#     gain      25000  50000  75000  100000  150000
#     hits/trk     46     68     77      80      89
# The sim par declares 150000, which would give 89 against a data median of 48.
#
# Four states run in parallel on 8 cores; the fifth follows. Each root job is single-threaded.
#
#   ./run_reco_batch_pt.sh [nEvents] [outDir]
set -uo pipefail
NEV=${1:-6000}
OUT=${2:-./data/reco}
PAR=${PAR:-ATTPC.a1975_C16_sim_g25000.par}
SIM=${SIM:-/mnt/f/a1975_C16_pt_sim}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../../../.." && pwd)
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include:${ROOT_INCLUDE_PATH:-}"; set -u
mkdir -p "$OUT"; cd "$HERE"

pids=()
for J in gs_s4001 ex1_s4001 ex2_s4001 ex3_s4001 ex4_s4001; do
  if [ -s "$OUT/${J}_reco.root" ]; then echo "  $J already reconstructed"; continue; fi
  [ -s "$SIM/${J}_sim.root" ] || { echo "  $J MISSING INPUT"; continue; }
  echo "[$(date +%H:%M:%S)] $J  ($NEV events, par $PAR)"
  ( root -b -l -q "run_reco_C16pt.C(\"$SIM/${J}_sim.root\",\"$OUT/${J}_reco.root\",\"$PAR\",20,20,8,$NEV)" \
      > "$OUT/${J}_reco.log" 2>&1
    grep -q "sim reco done" "$OUT/${J}_reco.log" && echo "  $J ok" || echo "  $J FAILED -- see $OUT/${J}_reco.log" ) &
  pids+=($!)
  # 4 at a time: 8 cores, and each job also drives I/O on /mnt/f
  while [ "$(jobs -rp | wc -l)" -ge 4 ]; do wait -n; done
done
wait
echo "[$(date +%H:%M:%S)] all states reconstructed in $OUT"
