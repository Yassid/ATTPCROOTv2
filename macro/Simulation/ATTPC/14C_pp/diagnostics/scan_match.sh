#!/usr/bin/env bash
# Scan Gain / CoefT / CoefL of the 14C(p,p') sim against the a1954 data reference.
#
# Reference (measured by match_diag on the surviving GATED data reco files,
# ~/a1954_C14_fit_300torr/in/, runs 0056 + 0061, same HDBSCAN config as the sim):
#     charge/mm  ~1070      hits/mm  ~1.72      charge/hit median ~280
#     sigma_xy^2 = 4.07 + 0.00085*z            sigma_z^2 = 0.20 + 0.00021*z
#
# Each point re-digitizes the SAME MC truth (data/attpcsim.root is fixed), so only
# Clusterize -> Pulse -> PSA -> clean -> PRA re-runs. Then gate_truth_C14.C selects the
# real recoil protons, because the data reference is IC+PID-gated and comparing gated
# data against an ungated sim (beam included) would be meaningless.
#
#   ./scan_match.sh [nEvents]      results -> scan_results.txt
set -eo pipefail

NEV=${1:-400}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
SIMDIR=$(cd "$HERE/.." && pwd)
REPO=$(cd "$SIMDIR/../../../.." && pwd)
PARDIR="$REPO/parameters"
BASEPAR="$PARDIR/ATTPC.a1954_C14_sim.par"
OUT="$HERE/scan"
RES="$HERE/scan_results.txt"
mkdir -p "$OUT"

set +u
source "$REPO/build/config.sh" >/dev/null 2>&1   # set -u + thisroot.sh kills the script
set -u

# tag:Gain:CoefT:CoefL
POINTS="
base:150000:0.0009:0.0009
T3:150000:0.0030:0.0002
T6:150000:0.0060:0.0002
T12:150000:0.0120:0.0002
T6L9:150000:0.0060:0.0009
T6g75:75000:0.0060:0.0002
"

echo "# scan started $(date +%H:%M:%S)  nEvents=$NEV" >> "$RES"
echo "# DATA REF: charge/mm 1070  hits/mm 1.72  q/hit(med) 280  sz2 = 0.20+0.00021z" >> "$RES"

for P in $POINTS; do
  TAG=${P%%:*}; REST=${P#*:}
  G=${REST%%:*}; REST=${REST#*:}
  CT=${REST%%:*}; CL=${REST##*:}

  PAR="ATTPC.a1954_C14_sim_$TAG.par"
  # NOTE: only ever REPLACE values on existing lines. A standalone comment line makes
  # FairParAsciiFileIo silently drop every parameter below it.
  sed -e "s/^Gain:Double_t.*/Gain:Double_t               $G   # scan $TAG/" \
      -e "s/^CoefT:Double_t.*/CoefT:Double_t              $CT   # scan $TAG/" \
      -e "s/^CoefL:Double_t.*/CoefL:Double_t              $CL   # scan $TAG/" \
      "$BASEPAR" > "$PARDIR/$PAR"
  grep -qE "^Gain:Double_t +$G" "$PARDIR/$PAR" || { echo "$TAG PAR_WRITE_FAILED" >> "$RES"; continue; }

  echo "[$(date +%H:%M:%S)] === $TAG  Gain=$G CoefT=$CT CoefL=$CL"
  RECO="$OUT/${TAG}_reco.root"
  ( cd "$SIMDIR" && root -b -q -l \
      "run_reco_C14.C(\"./data/attpcsim.root\",\"$RECO\",\"$PAR\",20,20,8,$NEV)" ) \
      > "$OUT/${TAG}_reco.log" 2>&1 || { echo "$TAG RECO_FAILED" >> "$RES"; continue; }
  [ -s "$RECO" ] || { echo "$TAG RECO_EMPTY" >> "$RES"; continue; }

  # truth-gate to the recoil protons (the data reference is gated)
  ( cd "$SIMDIR" && root -b -q -l \
      "gate_truth_C14.C(\"$RECO\",\"$OUT/\",\"${TAG}g\",0.6,4)" ) \
      > "$OUT/${TAG}_gate.log" 2>&1 || { echo "$TAG GATE_FAILED" >> "$RES"; continue; }
  GF="$OUT/${TAG}g_reco.root"
  [ -s "$GF" ] || { echo "$TAG GATE_EMPTY" >> "$RES"; continue; }

  { echo "### $TAG  Gain=$G CoefT=$CT CoefL=$CL"
    ( cd "$HERE" && ./match_diag "$GF" "$TAG" 1.30 -1 ) \
      | grep -E "tracks:|CHARGE|LENGTH|PAD-PLANE|DRIFT"
    echo "COMPLETED $TAG"; } >> "$RES"
  echo "[$(date +%H:%M:%S)] $TAG done"
done
echo "# scan finished $(date +%H:%M:%S)" >> "$RES"
