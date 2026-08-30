#!/usr/bin/env bash
# The 14C(d,p)15C field x pad-pitch campaign: 3 fields x 2 pitches x 2 levels.
#
#   ./run_C14dp_campaign.sh [-j N] [nEvents]
#
# 18 samples, 9 generations (both pitches of a field read the same sims, so a pitch
# difference is the pad plane and nothing else). Levels: the 1/2+ ground state and the 5/2+ at
# 0.740 MeV, which are the only two bound states of 15C.
#
# Wave 1 is the ground state across all six configurations, so the headline number exists early;
# wave 2 adds the 0.740, which is what turns a resolution into a separation.
set -uo pipefail
JOBS=5
if [ "${1:-}" = "-j" ]; then JOBS=${2:?-j needs a number}; shift 2; fi
NEV=${1:-8000}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOTDIR=${DP_ROOT:-/mnt/f/a1954_C14dp_hf}
mkdir -p "$ROOTDIR"
MASTER="$ROOTDIR/campaign.log"
say() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$MASTER"; }

STATES="gs:1 ex0740:2 ex3103:3"
FIELDS="2.85:0 4.0:1 7.0:2"
PADS="-1 2.0"

emit() {
   local want="$1"
   for spec in $STATES; do
      local st=${spec%%:*} off=${spec##*:}
      case " $want " in *" $st "*) ;; *) continue ;; esac
      for fspec in $FIELDS; do
         local b=${fspec%%:*} fi=${fspec##*:}
         local seed=$((8000 + 100 * fi + off))
         for pad in $PADS; do echo "$st:$b:$pad:$seed"; done
      done
   done
}
run_wave() {
   local label=$1 want=$2 list
   list=$(emit "$want")
   say "=== wave $label : $(echo "$list" | wc -l) samples, -j $JOBS, $NEV events each ==="
   printf "%s\n" $list | xargs -P "$JOBS" -I{} bash -c '
     IFS=: read -r st b pad seed <<< "{}"
     '"$HERE"'/accumulate_C14dp.sh "$st" "$b" "$pad" "$seed" '"$NEV"' 2>&1 | tail -2'
   say "=== wave $label finished ==="
}
say "########## (d,p) campaign start (-j $JOBS, $NEV events/sample) ##########"
run_wave 1 "gs ex0740"
run_wave 2 "ex3103"
say "########## (d,p) campaign done ##########"
for cfg in b285_attpc b285_2mm b400_attpc b400_2mm b700_attpc b700_2mm; do
   n=$(ls "$ROOTDIR/$cfg"/*.marker 2>/dev/null | wc -l)
   say "$(printf '%-14s %s/2' "$cfg" "$n")"
done
