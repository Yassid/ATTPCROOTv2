#!/usr/bin/env bash
# Reconstruct the SAME events twice, changing ONLY the pattern-recognition stage, so the
# clusterer can be blamed or cleared for the proton KE bias.
#
#   tc       = AtTrackFinderTC (triplclust)  -- what Ayyad et al., EPJ A 59:294 (2023) used
#   hdbscan  = AtTrackFinderHDBSCAN          -- what the current a1954 production uses
#
# Everything else is held fixed: same runs, same event range, same PSA (multifit, thr 20), same
# cleaner, same parameter file. Output goes to its own tree so nothing can overwrite production.
#
# NOTE the raw data moved from F: to G: (see reference_g_drive_a1954); the macros still default
# to /mnt/f/a1954_remerged, so the path is passed explicitly here.
#
#   ./pra_cmp_C14.sh "run_0058 run_0061" 12000

set -uo pipefail

RUNS="${1:-run_0058}"
NEV="${2:-12000}"
IN="/mnt/g/a1954_remerged/"
OUT="/home/yassid/a1954_C14_pracmp"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"

set +u
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
set -u

if [[ ! -d "$IN" ]]; then
   echo "FATAL: raw data dir $IN not found" >&2
   exit 1
fi

for pra in tc hdbscan; do
   mkdir -p "$OUT/$pra/logs"
done

for r in $RUNS; do
   for pra in tc hdbscan; do
      d="$OUT/$pra"
      # the marker is written only after the output file exists AND the macro exited 0, so a
      # killed job never reports success (see feedback_silent_degradation)
      if [[ -f "$d/${r}.marker" ]]; then
         echo "skip $r/$pra (already done)"
         continue
      fi
      echo "=== $r  PRA=$pra  nEv=$NEV ==="
      root -b -q -l "$HERE/pipeline/unpackReco_C14.C(\"$r\",$NEV,false,\"$d/\",\"$IN\",false,false,\"multifit\",0,20,\"$pra\")" \
         > "$d/logs/${r}_reco.log" 2>&1
      rc=$?
      if [[ $rc -eq 0 && -s "$d/${r}_reco.root" ]]; then
         touch "$d/${r}.marker"
         echo "  done $r/$pra"
      else
         echo "  FAILED $r/$pra (rc=$rc) -- see $d/logs/${r}_reco.log" >&2
      fi
   done
done
echo "ALL RECO DONE"
