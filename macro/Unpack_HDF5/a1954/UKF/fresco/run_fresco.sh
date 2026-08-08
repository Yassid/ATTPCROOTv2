#!/usr/bin/env bash
# Run FRESCO on one input from inputs/ and extract per-state dσ/dΩ curves to outputs/.
#
# Layout:
#   inputs/<stem>.nin       FRESCO input (version-controlled)
#   outputs/<stem>.out      full FRESCO text log
#   outputs/<stem>_dsdo.dat state 1 (elastic)   mb/sr vs deg
#   outputs/<stem>_dsdo_ex<N>.dat  state N>=2   mb/sr vs deg
#
# Usage: ./run_fresco.sh <input-stem>
#        ./run_fresco.sh p15C_el_13MeV
#        ./run_fresco.sh p15C_inel_13MeV_L2

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 <input-stem>" >&2
    exit 1
fi

STEM="${1%.nin}"
STEM="$(basename "$STEM")"

HERE="$(cd "$(dirname "$0")" && pwd)"
IN_DIR="$HERE/inputs"
OUT_DIR="$HERE/outputs"
IN="$IN_DIR/$STEM.nin"
OUT="$OUT_DIR/$STEM.out"

if [[ ! -f "$IN" ]]; then
    echo "error: input not found: $IN" >&2
    exit 1
fi
if ! command -v fresco >/dev/null; then
    echo "error: 'fresco' not in PATH (try ~/.local/bin)" >&2
    exit 1
fi

mkdir -p "$OUT_DIR"

# Run FRESCO in a scratch dir so fort.* don't pollute inputs/outputs.
SCRATCH="$(mktemp -d)"
trap 'rm -rf "$SCRATCH"' EXIT
cp "$IN" "$SCRATCH/"
( cd "$SCRATCH" && fresco < "$STEM.nin" > "$STEM.out" )

cp "$SCRATCH/$STEM.out" "$OUT"

# Split text output into per-state dat files.
awk -v outdir="$OUT_DIR" -v stem="$STEM" '
function flush() {
    if (state > 0 && np > 0) {
        if (state == 1) fn = outdir "/" stem "_dsdo.dat"
        else            fn = outdir "/" stem "_dsdo_ex" state ".dat"
        for (i = 0; i < np; i++) print lines[i] > fn
        close(fn)
        printf "  state %2d -> %s (%d angles)\n", state, fn, np
    }
    np = 0
}
/CROSS SECTIONS FOR OUTGOING.*state *# *[0-9]+/ {
    flush()
    match($0, /state *# *([0-9]+)/, m)
    state = m[1] + 0
    next
}
state > 0 && match($0, /^[[:space:]]*([0-9.]+)[[:space:]]*deg.*X-S =[[:space:]]*([0-9.Ee+-]+)/, m) {
    lines[np++] = m[1] " " m[2]
}
END { flush() }
' "$OUT"
