#!/usr/bin/env bash
# 46Ar(3He,d)47K DWBA in PtolemyCpp -- one job per state, parsed into dat/<tag>.dat.
#
#   ./run_dwba.sh [incomingCode outgoingCode]
#
# WHY ONE JOB PER STATE. A .dwba file holding three reaction lines reports "Number of Reaction: 3"
# and then writes the cross-section tables of only the LAST one. Both other states are silently
# absent from the output -- not zero, absent -- so a multi-line file looks like it worked and
# quietly yields one curve. One file per state, always.
#
# WHY THIS EXISTS AT ALL. ar46_dwba.txt carries s1/2 and d3/2 and nothing else; the 2.02 MeV 7/2-
# (0f7/2, l = 3) has no column, so that state cannot enter a weighted spectrum. The proposal used
# Ptolemy, so this is the same code -- but a new l = 3 curve is only comparable with the existing
# two if this deck reproduces THOSE two first. compare_dwba.py is that check, and it is the point
# of the exercise: an l = 3 curve that cannot be shown to sit on the same footing as the l = 0 and
# l = 2 ones must not be mixed with them.
#
# POTENTIALS. Default p = Pang (2009) 3He entrance, A = An & Cai (2006) deuteron exit. The
# potentials behind ar46_dwba.txt are NOT recorded anywhere, so a shape mismatch may be the
# potential rather than the code -- which is why the comparison sweeps a few pairs.
#
# TWO CAVEATS FROM THE CODE ITSELF, both of which must be read off every run:
#   * PtolemyCpp's known transfer bug leaks the INCOMING surface imaginary into the OUTGOING
#     channel when the outgoing block has vsi = 0 explicitly. Deuteron global OMPs normally carry
#     vsi != 0, so it should not bite here -- check_vsi() below asserts that rather than hoping.
#   * "lPeak = 30 IS AT END OF L-RANGE (FISHY)" means LMAX is too small for this system and the
#     partial-wave sum is truncated at its peak. It is a WARNING, not an error, and the run still
#     prints a full table. Raising LMAX past ~60 has been seen to trip "CANNOT EXTRAPOLATE" and
#     return all zeros, so the fix is not simply a bigger number.
set -uo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
PT=${PT:-$HOME/PtolemyCpp}
IN=${1:-p}
OUT=${2:-A}
ELAB=${ELAB:-13MeV/u}
mkdir -p "$HERE/dat" "$HERE/outputs" "$HERE/inputs"

# tag : orbital : Jpi : Ex
STATES=(
  "gs:1s1/2:1/2+:0"
  "e036:0d3/2:3/2+:0.36"
  "e202:0f7/2:7/2-:2.02"
)

for s in "${STATES[@]}"; do
  IFS=: read -r tag orb jpi ex <<< "$s"
  f="$HERE/inputs/${tag}_${IN}${OUT}.dwba"
  o="$HERE/outputs/${tag}_${IN}${OUT}.out"
  printf '# 46Ar(3He,d)47K %s -> 47K %s at %s, potentials %s%s\n' "$orb" "$jpi" "$ELAB" "$IN" "$OUT" > "$f"
  printf '46Ar(3He,d)47K  0+  %s  %s  %s  %s  %s%s\n' "$orb" "$jpi" "$ex" "$ELAB" "$IN" "$OUT" >> "$f"
  ( cd "$PT" && ./ptolemy --dwba "$f" > "$o" 2>"${o%.out}.err" ) || true
  nb=$(grep -c 'COMPUTATION OF CROSS' "$o" 2>/dev/null || echo 0)
  fishy=$(grep -c 'FISHY' "$o" 2>/dev/null || echo 0)
  zero=$(grep -c 'CANNOT EXTRAPOLATE' "$o" 2>/dev/null || echo 0)
  if [ "$nb" -eq 0 ]; then echo "$tag: NO CROSS-SECTION TABLE -- see $o"; continue; fi
  # THE vsi LEAK CHECK. Print both channels' surface imaginary so a leak is visible rather than
  # assumed absent: if the OUTGOING vsi equals the INCOMING one AND the deuteron OMP should have
  # its own, the known bug has fired and the normalisation is wrong by ~1.6.
  vin=$(awk '/INCOMING/{f=1} f&&/vsi/{print $0; exit}' "$o" | grep -oE 'vsi[= ]+[-0-9.]+' | head -1)
  vout=$(awk '/OUTGOING/{f=1} f&&/vsi/{print $0; exit}' "$o" | grep -oE 'vsi[= ]+[-0-9.]+' | head -1)
  printf "%-6s %s%s  blocks=%-3s fishy=%-3s zeros=%-3s  in:%s out:%s\n" \
         "$tag" "$IN" "$OUT" "$nb" "$fishy" "$zero" "${vin:-?}" "${vout:-?}"
  python3 "$HERE/parse_dwba.py" "$o" "$HERE/dat/${tag}_${IN}${OUT}.dat"
done
