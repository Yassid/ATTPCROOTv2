#!/bin/bash
# Step 2b -- TRACKING with ATTPCROOT's own pattern recognition (alternative to 02_track.py).
#
# Runs SampleConsensus::AtSampleConsensus, the method of Ayyad et al., NIM A 880 (2018) 166.
# Needs the ATTPCROOT environment, because this one IS the framework:
#     source <ATTPCROOT>/setup_fr19port.sh
#
#   ./02b_track_ransac.sh            all runs
#   ./02b_track_ransac.sh 3          first 3 runs, to try it
#
# Writes tracks/ransac_XXXX.csv, then converts them to the same pickle format 02_track.py
# produces, so steps 3-6 read either without changes.  Keep the two sets in different
# DATA_ROOTs if you want to compare them.
cd "$(dirname "$0")"
command -v root >/dev/null || { echo "ROOT not on PATH. source setup_fr19port.sh first."; exit 1; }
root -l -b -q 'lib/track_ransac.C' 2>&1 | grep -q usage || {
  echo "track_ransac.C will not load. Did you source setup_fr19port.sh (not just thisroot.sh)?"; exit 1; }
python3 - <<'PY' > /tmp/dec_paths.txt
import sys; sys.path.insert(0,'.')
import config; print(config.ROOT_DIR); print(config.TRACKS_DIR)
PY
ROOT_DIR=$(sed -n 1p /tmp/dec_paths.txt); TRACKS=$(sed -n 2p /tmp/dec_paths.txt)
mkdir -p "$TRACKS"
ITER=${ITER:-1000}; MINHITS=${MINHITS:-10}; DIST=${DIST:-6}
echo "  settings: iterations=$ITER minHits=$MINHITS distanceThreshold=$DIST mm"
N=${1:-0}; i=0
for f in "$ROOT_DIR"/alpha_run_*_hits.root; do
  [ -e "$f" ] || { echo "no ROOT files in $ROOT_DIR"; exit 1; }
  run=$(basename "$f" | sed 's/alpha_run_\([0-9]*\)_hits.root/\1/')
  csv="$TRACKS/ransac_$run.csv"
  [ -s "$csv" ] || root -l -b -q "lib/track_ransac.C(\"$f\",\"$csv\",-1,$ITER,$MINHITS,$DIST)" 2>&1 | grep -E "tracks ->|ERROR"
  python3 lib/ransac_to_pickle.py "$csv" "$TRACKS/tracks_$run.pkl"
  i=$((i+1)); [ "$N" -gt 0 ] && [ "$i" -ge "$N" ] && break
done
echo "  tracks in $TRACKS"
