#!/usr/bin/env bash
# The field x pad-pitch matrix for 46Ar(3He,d)47K, run end to end.
#
#   ./run_3Hed_campaign.sh [-j N]
#
#   config              B [T]   pads      sims from              output
#   A (already done)     2.85   AT-TPC    /mnt/f/ar46_3hed       /mnt/f/ar46_3hed
#   B                    3.80   AT-TPC    generated here         /mnt/f/ar46_3hed_B38
#   C                    2.85   2 mm      A's sims (reused)      /mnt/f/ar46_3hed_2mm
#   D                    3.80   2 mm      B's sims (reused)      /mnt/f/ar46_3hed_B38_2mm
#
# ONLY ONE NEW GENERATION SET IS PRODUCED. Transport depends on the field but not on the pad
# plane, so C reads the 2.85 T sims that already exist and D reads the ones B writes. That saves
# twelve generations and, more importantly, means each field's two pad configurations see the
# SAME events -- so a difference between them is the pad plane and nothing else.
#
# ORDER MATTERS: D consumes B's sims, so B runs first. The configurations run sequentially rather
# than all at once; -j applies within a configuration.
#
# COST, from the measured 0.42 s/entry (the same for both pad pitches -- 2 mm is not slower):
# ~90 min per 12000-entry sample, six samples per configuration, three configurations at -j 3
# is roughly 9 hours.
set -uo pipefail
JOBS=3
if [ "${1:-}" = "-j" ]; then JOBS=${2:?-j needs a number}; shift 2; fi
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MASTER=/mnt/f/ar46_3hed_campaign.log

say() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$MASTER"; }

run_cfg() { # name  BT  PAD  SIMDIR  OUT
   local name=$1 bt=$2 pad=$3 simdir=$4 out=$5
   say "=== config $name : B = $bt T, pads = $pad mm, sims $simdir -> $out ==="
   local t0=$SECONDS
   OUT="$out" BT="$bt" PAD="$pad" SIMDIR="$simdir" "$DIR/run_3Hed_accumulation.sh" -j "$JOBS"
   local n=$(ls "$out"/*.marker 2>/dev/null | wc -l)
   say "=== config $name finished: $n/6 samples, $(( (SECONDS-t0)/60 )) min ==="
}

say "########## campaign start (-j $JOBS) ##########"
run_cfg B 3.8  -1  /mnt/f/ar46_3hed_B38      /mnt/f/ar46_3hed_B38
run_cfg C 2.85 2.0 /mnt/f/ar46_3hed          /mnt/f/ar46_3hed_2mm
run_cfg D 3.8  2.0 /mnt/f/ar46_3hed_B38      /mnt/f/ar46_3hed_B38_2mm
say "########## campaign done ##########"
for d in /mnt/f/ar46_3hed /mnt/f/ar46_3hed_B38 /mnt/f/ar46_3hed_2mm /mnt/f/ar46_3hed_B38_2mm; do
   say "$(basename "$d"): $(ls "$d"/*.marker 2>/dev/null | wc -l)/6 complete"
done
