#!/bin/bash
# Build symlinks with contiguous file<N>_ indices (AtGRAWUnpacker maps files to decoders by
# a "%i" pattern) and a run list ordered chunk-0-first (ls sorts ".1.graw" before ".graw",
# which desynchronises the per-CoBo event counters and silently kills PSA).
#
# Writes the list to BOTH the scratch area and the canonical location the unpack drivers
# read, $VMCWORKDIR/runfiles/NSCL/Dec2014_alphas/alpha_<run>.txt. Those drivers used to
# look for a file this script never produced (different directory AND an "alpha_" prefix
# this script omitted), so a fresh checkout reported "[skip] : no runfile" for every run.
#
# Override the two roots if your data or checkout live elsewhere:
#   ATTPC_RAW   dir holding run_XXXX/CoBo*.graw   (default: the Cris_OneD drive)
#   ATTPC_LINKS scratch dir for the symlink trees
#
# Prints: <runfile path> <number of CoBos>
run=$1
RAW=${ATTPC_RAW:-/media/yassid/Cris_OneD/Dec2014_alphas}
LINKS=${ATTPC_LINKS:-/home/yassid/dec2014_links}
WORK=${VMCWORKDIR:-/home/yassid/fair_install/ATTPCROOTv2_fr19port}
D=$RAW/$run
L=$LINKS/$run
[ -d "$D" ] || { echo "NO_DIR"; exit 1; }
rm -rf "$L"; mkdir -p "$L"
cobos=$(ls "$D"/CoBo*.graw 2>/dev/null | sed 's|.*/CoBo\([0-9]*\)_.*|\1|' | sort -n -u)
i=0
for c in $cobos; do
  found=0
  for f in $(ls "$D"/CoBo${c}_*.graw 2>/dev/null | perl -pe 's|(.*?(?:\.([0-9]+))?\.graw)|($2//0)."\t".$1|e' | sort -n -k1,1 | cut -f2); do
    [ -s "$f" ] || continue
    ch=$(basename "$f" | grep -oE '\.[0-9]+\.graw$' | grep -oE '[0-9]+'); ch=${ch:-0}
    ln -s "$f" "$(printf "%s/file%d_c%02d_%s" "$L" "$i" "$ch" "$(basename "$f")")"
    found=1
  done
  [ $found -eq 1 ] && i=$((i+1)) || true
done
[ $i -eq 0 ] && { echo "NO_DATA"; exit 1; }
CANON="$WORK/runfiles/NSCL/Dec2014_alphas/alpha_${run}.txt"
mkdir -p "$(dirname "$CANON")"
ls "$L"/* > "$CANON"
cp "$CANON" "$LINKS/$run.txt"    # legacy location, kept so older drivers still work
echo "$CANON $i"
