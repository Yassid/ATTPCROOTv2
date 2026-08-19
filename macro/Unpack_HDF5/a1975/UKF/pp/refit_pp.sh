#!/usr/bin/env bash
# Refit the a1975 (p,p) runs with a different proton PID gate.
#
# Only the FIT stage is rerun. The reco files (<run>_reco.root, 165 GB of them) are untouched and
# AtPIDTask writes the PID landscape for every pattern track regardless of the gate, so the gate
# changes exactly one thing: which tracks genfit is asked to fit. That is the whole reason the
# existing production cannot simply be re-cut offline -- the tracks the old gate rejected were
# never fitted, so there is no Ex for them to recover.
#
# Output goes to <run>_genfitter<suffix>.root, so the existing _pp production stays where it is and
# the two can be compared run by run.
#
#   ./refit_pp.sh <gateJson> <suffix> [nParallel] [outDir]
#   ./refit_pp.sh pid/proton_sim.json _ppsim 4
#
# The gate path is relative to macro/Unpack_HDF5/a1975/UKF, the directory the fit macro resolves
# from. A marker file is written only after the output exists AND is non-empty -- a written file is
# not a finished job.
set -eo pipefail
GATE=${1:?need a gate json, e.g. pid/proton_sim.json}
SUF=${2:?need an output suffix, e.g. _ppsim}
NPAR=${3:-4}
OUT=${4:-/mnt/f/a1975/reco/}
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
UKF=$REPO/macro/Unpack_HDF5/a1975/UKF
MARK=/mnt/f/a1975_pp_refit${SUF}
mkdir -p "$MARK"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
cd "$UKF"
[ -f "$GATE" ] || { echo "gate $GATE not found (paths are relative to $UKF)"; exit 1; }

one() {
  R=$1
  [ -f "$MARK/$R.marker" ] && return 0
  root -b -q -l "pipeline/fitGenfitter_a1975.C(\"$R\",-1,\"/mnt/f/a1975/reco/\",\"$SUF\",\"$OUT\",-2.85,2,5,\"proton_H2_catima.txt\",\"$GATE\",4.0,kFALSE,2212,1.00727647,1,kFALSE,\"ATTPC_H300torr_RT\")" \
       > "$MARK/$R.log" 2>&1
  # the gate silently disables itself if the path does not resolve -- that would produce a full
  # ungated refit that looks like success, so check for it explicitly
  grep -q "gating will be disabled" "$MARK/$R.log" && { echo "  $R GATE_NOT_FOUND"; return 1; }
  [ -s "${OUT}${R}_genfitter${SUF}.root" ] || { echo "  $R FIT_FAILED"; return 1; }
  # A NON-EMPTY file is not a finished job either. run_0158 hit "Cannot allocate memory" reading
  # its 931 MB reco file, so TChain never found cbmsim, nothing was processed, and the 18 kB file
  # it left behind -- an empty cbmsim -- sailed through -s and got a COMPLETED marker. It was only
  # caught eight days later, in the cache. Count the entries; the run took 3.8 s instead of minutes.
  grep -q "Cannot allocate memory" "$MARK/$R.log" && { echo "  $R OOM_READING_RECO -- rerun it alone"; return 1; }
  NENT=$(root -b -q -l -e "TFile f(\"${OUT}${R}_genfitter${SUF}.root\"); TTree *t=(TTree*)f.Get(\"cbmsim\"); printf(\"NENT=%lld\\n\", t ? t->GetEntries() : -1);" 2>/dev/null \
         | sed -n 's/.*NENT=\(-\?[0-9]*\).*/\1/p' | tail -1)
  [ "${NENT:-0}" -gt 0 ] 2>/dev/null || { echo "  $R FIT_EMPTY -- cbmsim has ${NENT:-no} entries"; return 1; }
  echo COMPLETED > "$MARK/$R.marker"
  echo "[$(date +%H:%M:%S)] $R done"
}
export -f one; export MARK SUF OUT GATE

RUNS=$(ls /mnt/f/a1975/reco/*_reco.root | sed 's|.*/||; s|_reco.root||' | sort)
echo "[$(date +%H:%M:%S)] refitting $(echo "$RUNS" | wc -l) runs, gate $GATE, suffix $SUF, $NPAR at a time"
echo "$RUNS" | xargs -P "$NPAR" -I{} bash -c 'one {}' || true
echo "[$(date +%H:%M:%S)] REFITDONE  $(ls "$MARK"/*.marker 2>/dev/null | wc -l) / $(echo "$RUNS" | wc -l) completed"
