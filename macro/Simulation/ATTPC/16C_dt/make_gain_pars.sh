#!/usr/bin/env bash
# Write Gain variants of the a1975 D2 production par, for the (d,t) simulation gain scan.
#
# WHY A SCAN IS NEEDED AT ALL. Gain decides how many electrons reach a pad, so it decides whether
# a pad crosses the PSA threshold, so it decides how many hits a track keeps -- which is what an
# acceptance measures. The production par declares Gain = 10000 while the (p,d) SIMULATION par
# declares 150000, and that one carries its own g400000/g1000000/g2500000/g5000000 variants, i.e.
# it was tuned rather than guessed. Using the data par's gain unexamined would be assuming that a
# number chosen for reconstructing real electronics also describes simulated charge collection.
#
# WHAT IS AND IS NOT VARIED. Only the Gain line changes. Drift velocity, TBEntrance, ZPadPlane,
# density, pressure and every geometry field stay exactly as the production par has them, because
# those must match the data or the simulation's z scale is not the data's z scale. The diff
# against the source par is printed for each variant so that claim is checkable, not asserted.
#
# HOW TO CHOOSE THE WINNER: not by which gain gives the highest acceptance -- that is monotonic
# and would always pick the largest. Compare HITS PER TRACK and PADS PER TRACK against the DATA
# (reco_d2_dv1104). The (p,d) chain is already known to diverge this way, 54.4 clusters/track in
# simulation against 34.5 in data, so this is a measured failure mode.
#
#   ./make_gain_pars.sh [gains...]
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../../../.." && pwd)
SRC="$REPO/parameters/ATTPC.a1975_deuterium_dv1104.par"
GAINS=${*:-"50000 150000 400000 1000000"}

[ -s "$SRC" ] || { echo "ERROR: missing $SRC"; exit 1; }
grep -qE '^Gain:' "$SRC" || { echo "ERROR: no Gain line in $SRC -- refusing to guess"; exit 1; }
echo "source: $SRC  (Gain = $(grep -E '^Gain:' "$SRC" | awk '{print $2}'))"

for g in $GAINS; do
  OUT="$REPO/parameters/ATTPC.a1975_deuterium_dv1104_g${g}.par"
  sed -E "s|^(Gain:[A-Za-z_]*[[:space:]]+)[0-9.eE+-]+|\1${g}|" "$SRC" > "$OUT"
  # Verify the substitution took AND that it is the ONLY difference. A sed that silently matched
  # nothing leaves a file that looks fine and quietly reruns the production gain.
  got=$(grep -E '^Gain:' "$OUT" | awk '{print $2}')
  ndiff=$(diff "$SRC" "$OUT" | grep -c '^[<>]')
  if [ "$got" != "$g" ]; then
    echo "  FAILED  g=$g: par still says $got"; rm -f "$OUT"; continue
  fi
  if [ "$ndiff" -ne 2 ]; then
    echo "  FAILED  g=$g: $ndiff changed lines, expected exactly 2 (the Gain line, before and after)"
    rm -f "$OUT"; continue
  fi
  echo "  ok  Gain=$g  -> $(basename "$OUT")"
done
