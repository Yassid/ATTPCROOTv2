#!/usr/bin/env bash
# Extra statistics for the BACKWARD region of 14C(d,p)15C, generated over theta_cm 2-60 deg only.
#
#   ./run_C14dp_backward.sh [-j N] [nEvents]
#
# The uniform 2-178 campaign puts only ~180 reconstructed protons per configuration into
# theta_cm 8-30 -- the transfer peak, and the slice the whole channel was simulated for. Restricting
# the generation to 2-60 deg concentrates every event there instead. The acceptance is a per-bin
# ratio and the resolution is measured per slice, so restricting the range costs nothing in
# validity; it only stops paying for events at angles that are not the question.
#
# Separate DP_ROOT so these never mix with the uniform samples in a glob.
set -uo pipefail
JOBS=5
if [ "${1:-}" = "-j" ]; then JOBS=${2:?-j needs a number}; shift 2; fi
NEV=${1:-8000}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
export DP_ROOT=${DP_ROOT:-/mnt/f/a1954_C14dp_back}
export CM_LO=2.0 CM_HI=60.0
mkdir -p "$DP_ROOT"
MASTER="$DP_ROOT/campaign.log"
say() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$MASTER"; }
say "########## (d,p) BACKWARD campaign, theta_cm $CM_LO-$CM_HI, -j $JOBS, $NEV events ##########"
LIST=""
for spec in gs:1 ex0740:2; do
   st=${spec%%:*}; off=${spec##*:}
   for fspec in 2.85:0 4.0:1 7.0:2; do
      b=${fspec%%:*}; fi=${fspec##*:}
      seed=$((8500 + 100 * fi + off))
      for pad in -1 2.0; do LIST="$LIST $st:$b:$pad:$seed"; done
   done
done
printf "%s\n" $LIST | xargs -P "$JOBS" -I{} bash -c '
  IFS=: read -r st b pad seed <<< "{}"
  '"$HERE"'/accumulate_C14dp.sh "$st" "$b" "$pad" "$seed" '"$NEV"' 2>&1 | tail -2'
say "########## backward campaign done ##########"
for cfg in b285_attpc b285_2mm b400_attpc b400_2mm b700_attpc b700_2mm; do
   say "$(printf '%-14s %s/2' "$cfg" "$(ls "$DP_ROOT/$cfg"/*.marker 2>/dev/null | wc -l)")"
done
