#!/usr/bin/env bash
# The 10Be(t,p)12Be field x pad-pitch campaign: 3 fields x 2 pitches x 4 levels = 24 samples.
#
#   ./run_Be10tp_campaign.sh [-j N] [nEvents]
#
# Levels: the four BOUND states of 12Be -- 0+ g.s., 2+ 2.109, 0+_2 2.251, 1- 2.715 (S_n = 3.171).
# All four are simulated at FULL statistics; the 5x weaker population of the 0+_2 is applied when
# the four are summed into one spectrum (tp_spectrum_Be10.C), NOT by generating fewer events. See
# the header of accumulate_Be10tp.sh for why.
#
# ONE xargs STREAM, not waves. The (d,p) campaign used wave barriers and paid for them in idle
# cores at each boundary; here the whole 24-sample list is fed to a single -P N pool so nothing
# idles. The list is ordered FIELD-MAJOR and, within a field, level-major with the two pad planes
# adjacent -- so (a) the 2.85 T baseline, which is the real detector, finishes first and the
# headline spectrum exists early, and (b) the two pitches of a level sit next to each other, so the
# second one finds the shared generation already done (or waits a couple of minutes on its lock)
# instead of duplicating it.
set -uo pipefail
JOBS=6
if [ "${1:-}" = "-j" ]; then JOBS=${2:?-j needs a number}; shift 2; fi
NEV=${1:-16000}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOTDIR=${TP_ROOT:-/mnt/f/Be10_tp}
mkdir -p "$ROOTDIR"
MASTER="$ROOTDIR/campaign.log"
say() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$MASTER"; }

STATES="gs:1 ex2109:2 ex2251:3 ex2715:4"
FIELDS="2.85:0 4.0:1 7.0:2"
PADS="-1 2.0"

emit() {
   for fspec in $FIELDS; do
      local b=${fspec%%:*} fi=${fspec##*:}
      for spec in $STATES; do
         local st=${spec%%:*} off=${spec##*:}
         local seed=$((9000 + 100 * fi + off))
         for pad in $PADS; do echo "$st:$b:$pad:$seed"; done
      done
   done
}

LIST=$(emit)
say "########## (t,p) campaign start (-j $JOBS, $NEV events/sample, $(echo "$LIST" | wc -l) samples) ##########"
printf "%s\n" $LIST | xargs -P "$JOBS" -I{} bash -c '
  IFS=: read -r st b pad seed <<< "{}"
  '"$HERE"'/accumulate_Be10tp.sh "$st" "$b" "$pad" "$seed" '"$NEV"' 2>&1 | tail -2'
say "########## (t,p) campaign done ##########"
for cfg in b285_attpc b285_2mm b400_attpc b400_2mm b700_attpc b700_2mm; do
   n=$(ls "$ROOTDIR/$cfg"/*.marker 2>/dev/null | wc -l)
   say "$(printf '%-14s %s/4' "$cfg" "$n")"
done
