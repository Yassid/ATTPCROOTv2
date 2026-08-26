#!/usr/bin/env bash
# Run every deck in inputs/ through PtolemyCpp and parse the angular distributions into dat/.
#
# PtolemyCpp (https://github.com/goluckyryan/PtolemyCpp) is a C++ reimplementation of the Argonne
# Ptolemy code. It needs nothing but g++ with C++17:
#     git clone https://github.com/goluckyryan/PtolemyCpp.git && cd PtolemyCpp && make
# Point PTOLEMY at the resulting binary, or keep it at the default below.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PTOLEMY="${PTOLEMY:-/home/yassid/PtolemyCpp/ptolemy}"
[ -x "$PTOLEMY" ] || { echo "no ptolemy binary at $PTOLEMY -- set PTOLEMY=/path/to/ptolemy"; exit 1; }
mkdir -p "$HERE/outputs" "$HERE/dat"
for f in "$HERE"/inputs/*; do
  b=$(basename "$f"); b=${b%.*}
  "$PTOLEMY" < "$f" > "$HERE/outputs/$b.out" 2>&1 || echo "  $b: nonzero exit"
  printf '  %-28s %s\n' "$b" "$(grep -c 'CANNOT EXTRAPOLATE' "$HERE/outputs/$b.out" || true) extrap errors"
done
python3 "$HERE/parse_ptolemy.py" "$HERE"
