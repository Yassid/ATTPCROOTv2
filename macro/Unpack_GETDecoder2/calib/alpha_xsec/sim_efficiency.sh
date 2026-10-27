#!/bin/bash
# Measure the tracking efficiency of a simulated sample against its own MC truth.
#
#     ./sim_efficiency.sh <real.root> <truth.csv>
#     TRACKER=ransac ./sim_efficiency.sh <real.root> <truth.csv>
#
# THE TRACKER IS THE KNOB.  Everything else in the chain -- the simulation, the
# digitisation, the merge, the pile-up, the reconstruction and the cuts -- is the same code
# the real data go through, so whatever you change about the tracking is the only thing the
# comparison is sensitive to.
#
#     TRACKER=hdbscan   (default)  02_track.py, pure Python, no framework needed
#     TRACKER=ransac               AtSampleConsensus from ATTPCROOT
#
# To try your OWN finder, write something that produces the same per-event pickle -- a dict
# with tb/x/y/q arrays and three boolean masks beam/a1/a2 plus the event number -- and point
# TRACKER at it.  See lib/ransac_to_pickle.py for the smallest possible example.
set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
REAL=${1:?usage: sim_efficiency.sh <real.root> <truth.csv>}
TRUTH=${2:?usage: sim_efficiency.sh <real.root> <truth.csv>}
TRACKER=${TRACKER:-hdbscan}
WORK=${WORK:-$HERE/out/simwork}
mkdir -p "$WORK/csv" "$WORK/tracks"

command -v root >/dev/null || { echo "ROOT is not on the PATH"; exit 1; }
echo "[$(date +%T)] exporting hits from $(basename "$REAL")"
root -l -b -q "lib/dump_hits.C(\"$REAL\",\"$WORK/csv/hits_9999.csv\")" 2>&1 | grep -E "hits ->|ERROR"

case "$TRACKER" in
  hdbscan)
    echo "[$(date +%T)] tracking with HDBSCAN [position, orientation]"
    DEC2014_DATA="$WORK" python3 02_track.py 2>&1 | tail -1
    PKL=$(ls "$WORK"/tracks/tracks_*.pkl | head -1)
    ;;
  ransac)
    echo "[$(date +%T)] tracking with AtSampleConsensus (needs the ATTPCROOT environment)"
    root -l -b -q "lib/track_ransac.C(\"$REAL\",\"$WORK/ransac.csv\")" 2>&1 | tail -1
    python3 lib/ransac_to_pickle.py "$WORK/ransac.csv" "$WORK/tracks/tracks_9999.pkl"
    PKL="$WORK/tracks/tracks_9999.pkl"
    ;;
  *)
    echo "unknown TRACKER=$TRACKER (expected hdbscan or ransac)"; exit 1 ;;
esac

echo "[$(date +%T)] scoring against MC truth"
python3 lib/score_tracking.py "$TRUTH" "$PKL" "$HERE/out/sim_efficiency.png" "$TRACKER"
echo "[$(date +%T)] done -- out/sim_efficiency.png and out/efficiency.npz"
echo "   out/efficiency.npz is what 07_corrected_xsec.py divides the cross section by."
