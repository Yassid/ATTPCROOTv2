#!/usr/bin/env bash
# The 14C(p,p') field x pad-pitch campaign: 3 fields x 2 pitches x 5 levels.
#
#   ./run_C14_hf_campaign.sh [-j N] [nEvents]
#
#   config        B [T]   pads       sims
#   b285_attpc     2.85   AT-TPC     sims_b285   <- the a1954 running conditions, the anchor
#   b285_2mm       2.85   2 mm       sims_b285
#   b400_attpc     4.00   AT-TPC     sims_b400
#   b400_2mm       4.00   2 mm       sims_b400
#   b700_attpc     7.00   AT-TPC     sims_b700
#   b700_2mm       7.00   2 mm       sims_b700
#
# Both pitches of a field READ THE SAME SIMS, so a pitch difference is the pad plane and nothing
# else; the fields necessarily have their own transport. Fifteen generations, thirty
# reconstructions.
#
# LEVELS: the elastic channel plus the four states the a1954 analysis actually has to separate --
# 6.094 (1-), 6.728 (3-), 7.012 (2+, the B(E2) carrier) and 8.317 (2+, isolated). The whole point
# of the matrix is whether a finer pitch or a higher field pulls 6.728 and 7.012 apart.
#
# RUN ORDER IS DELIBERATE. Wave 1 does gs and 6.094 across all six configurations, so the
# headline resolution comparison exists before the campaign is half over; wave 2 fills in the
# three remaining levels. Everything is resumable, so a wave that is interrupted costs only the
# sample it was in the middle of.
set -uo pipefail
JOBS=4
if [ "${1:-}" = "-j" ]; then JOBS=${2:?-j needs a number}; shift 2; fi
NEV=${1:-8000}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOTDIR=${HF_ROOT:-/mnt/f/a1954_C14_hf}
mkdir -p "$ROOTDIR"
MASTER="$ROOTDIR/campaign.log"

say() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$MASTER"; }

# state:seed-offset. The seed is 7000 + 100*field-index + offset, so no two samples anywhere in
# the matrix share a random sequence.
STATES="gs:1 ex6094:2 ex6728:3 ex7012:4 ex8317:5"
FIELDS="2.85:0 4.0:1 7.0:2"
PADS="-1 2.0"

emit_wave() { # $1 = space-separated state tags to include
   local want="$1" st off b fi pad seed s
   for spec in $STATES; do
      st=${spec%%:*}; off=${spec##*:}
      case " $want " in *" $st "*) ;; *) continue ;; esac
      for fspec in $FIELDS; do
         b=${fspec%%:*}; fi=${fspec##*:}
         seed=$((7000 + 100 * fi + off))
         for pad in $PADS; do
            echo "$st:$b:$pad:$seed"
         done
      done
   done
}

run_wave() { # $1 = label, $2 = states
   local label=$1 want=$2
   local list; list=$(emit_wave "$want")
   say "=== wave $label : $(echo "$list" | wc -l) samples, -j $JOBS, $NEV events each ==="
   printf "%s\n" $list | xargs -P "$JOBS" -I{} bash -c '
     IFS=: read -r st b pad seed <<< "{}"
     '"$HERE"'/accumulate_C14_hf.sh "$st" "$b" "$pad" "$seed" '"$NEV"' 2>&1 | tail -2'
   say "=== wave $label finished ==="
}

say "########## campaign start (-j $JOBS, $NEV events/sample) ##########"
run_wave 1 "gs ex6094"
say "--- wave 1 status ---"; "$HERE"/hf_status.sh | tee -a "$MASTER"
run_wave 2 "ex6728 ex7012 ex8317"
say "########## campaign done ##########"
"$HERE"/hf_status.sh | tee -a "$MASTER"
