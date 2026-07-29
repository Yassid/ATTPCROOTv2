#!/bin/bash
# Move the existing UKF fits to the external drive, then run the full three-fitter
# pass there. /home cannot hold the fit output (~140 GB for 39 runs x 3 fitters, see
# a2091 notes); the Seagate has 1.7 TB. Migrating rather than recomputing also saves
# ~3 h of CPU, and it makes fitall skip the 28 UKF fits it would otherwise redo.
set -o pipefail
HERE=/home/yassid/fair_install/ATTPCROOTv2/macro/Unpack_HDF5/a2091/UKF
SRC=/home/yassid/a2091_C15_fit
DST="/media/yassid/Seagate Hub/ATTPC/a2091_C15_fit"
LOG=/home/yassid/a2091_C15_reco/logs
say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*"; }

mkdir -p "$DST"
n=$(ls "$SRC"/*_ukf.root 2>/dev/null | wc -l)
say "migrating $n UKF fits -> external"
rsync -a --info=progress2 "$SRC"/*_ukf.root "$DST"/ || { say "RSYNC FAILED - keeping originals"; exit 1; }

# verify every file arrived at identical size BEFORE deleting anything
bad=0
for f in "$SRC"/*_ukf.root; do
  b=$(basename "$f")
  s1=$(stat -c%s "$f"); s2=$(stat -c%s "$DST/$b" 2>/dev/null || echo 0)
  [ "$s1" = "$s2" ] || { say "SIZE MISMATCH $b ($s1 vs $s2)"; bad=1; }
done
[ "$bad" -ne 0 ] && { say "verification failed - originals kept, aborting"; exit 1; }
say "verified $n files; removing local copies"
rm -f "$SRC"/*_ukf.root
say "/home now: $(df -BG --output=avail /home | tail -1 | tr -dc '0-9')GB free"

# Full pass over every reco'd run (auto-discovered), output on the external.
# UKF is skipped where it already exists; GENFIT noMat + matFX are what we need.
say "launching full three-fitter pass -> external"
FIT="$DST" "$HERE/fitall_C15.sh" "" 6 >> "$LOG/fitall_seagate.log" 2>&1
say "DONE: ukf=$(ls "$DST"/*_ukf.root 2>/dev/null | wc -l) nomat=$(ls "$DST"/*_genfit_nomat.root 2>/dev/null | wc -l) matFX=$(ls "$DST"/*_genfit_mat.root 2>/dev/null | wc -l)"
