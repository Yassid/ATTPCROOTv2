#!/usr/bin/env bash
# Rebuild the C15d PID explorer from the cached planes and open it.
#
#   ./open_pid_explorer.sh                 # D2 runs 17-103, gain matched
#   ./open_pid_explorer.sh raw             # same plane, no gain table (for comparison)
#   ./open_pid_explorer.sh d2 17 103       # explicit run range
#   ./open_pid_explorer.sh h2              # the EXCLUDED hydrogen runs, for inspection only
#
# The page is one self-contained HTML file with the data baked in: no server, no ROOT and no X11
# needed to VIEW it. X11 is only used here to launch the browser; the file can equally be copied
# somewhere and opened by hand.
#
# The `h2` mode exists so the hydrogen runs can be looked at, not analysed. They are a different
# target (dE/dx halves across run 103/106) and must never be pooled with the D2 set -- the page
# title says so when you select it.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"

MODE="${1:-d2}"
case "$MODE" in
   d2)  RMIN="${2:-17}";  RMAX="${3:-103}"; GAIN="gainmatch_C15d.csv"; TAG="15C + d  --  PID plane (D2)";;
   raw) RMIN="${2:-17}";  RMAX="${3:-103}"; GAIN="";                   TAG="15C + d  --  PID plane (D2, RAW)";;
   h2)  RMIN="${2:-106}"; RMAX="${3:-133}"; GAIN="";                   TAG="a1975 H2 runs -- NOT the 15C+d analysis";;
   *)   echo "usage: $0 [d2|raw|h2] [runMin] [runMax]" >&2; exit 1;;
esac

OUT="$HOME/C15d_pid_explorer${MODE:+_$MODE}.html"
[[ "$MODE" == d2 ]] && OUT="$HOME/C15d_pid_explorer.html"

# config.sh reads unset vars: source it with `set +u` or it kills the script silently
set +u
# shellcheck disable=SC1091
source "$REPO/build/config.sh" >/dev/null 2>&1
set -u

echo "building $OUT  (mode $MODE, runs $RMIN-$RMAX)"
root -b -q "$HERE/make_pid_explorer_html.C(\"/home/yassid/C15d_reco/\",\"$OUT\",\"$GAIN\",$RMIN,$RMAX)" \
   2>&1 | grep -vE '^Info|^Processing|^\s*\||^\s*-----|Welcome'

[[ -s "$OUT" ]] || { echo "ERROR: $OUT was not written" >&2; exit 1; }

if [[ -n "${DISPLAY:-}${WAYLAND_DISPLAY:-}" ]]; then
   echo "opening in the default browser ..."
   # Detach fully: the browser must outlive this script, and its chatter must not land in the
   # caller's terminal.
   setsid xdg-open "$OUT" >/dev/null 2>&1 < /dev/null &
   disown || true
else
   echo "no display -- open $OUT yourself, or copy it somewhere with a browser."
fi
