#!/usr/bin/env bash
# Open the C15d gate drawer.
#
#   ./pid/open_gate_gui.sh                          # proton gate, no IC window, no locus
#   ./pid/open_gate_gui.sh deuteron_C15d            # name the output gate
#   ./pid/open_gate_gui.sh proton_C15d 1000 1300    # with an IC window
#
# Run it from the workspace root (macro/Unpack_HDF5/C15d), or from anywhere -- it cd's there
# itself, because the macro resolves pid/points_C15d.root and pid/plots/ relative to the cwd.
#
# It is INTERACTIVE: `root -l`, never -b. It needs a display.

set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"     # .../C15d
REPO="$(cd "$HERE/../../.." && pwd)"

NAME="${1:-proton_C15d}"
ICLO="${2:--1}"
ICHI="${3:--1}"
EBEAM="${4:-0}"          # 0 disables the kinematic locus; pass a real lab energy once calibrated

# config.sh reads unset variables: source it under `set +u` or it exits the shell silently.
set +u
# shellcheck disable=SC1091
source "$REPO/build/config.sh" >/dev/null 2>&1
set -u

command -v root >/dev/null || { echo "ERROR: root not on PATH after sourcing config.sh" >&2; exit 1; }
[[ -n "${DISPLAY:-}${WAYLAND_DISPLAY:-}" ]] || { echo "ERROR: no display -- this GUI needs one" >&2; exit 1; }

cd "$HERE"
[[ -s pid/points_C15d.root ]] || {
   echo "ERROR: pid/points_C15d.root missing. Build it first:" >&2
   echo "   root -b -q 'pid/make_points_C15d.C()'" >&2
   exit 1
}
mkdir -p pid/plots

echo "opening the gate drawer: gate=$NAME  IC=[$ICLO,$ICHI]  Ebeam=$EBEAM"
# ★ ROOT's interactive prompt QUITS ON STDIN EOF. Launched detached with stdin from /dev/null the
# GUI appears and vanishes in the same second, which looks exactly like a crash. Keep stdin open.
[[ -t 0 ]] || exec 0< <(while :; do sleep 3600; done)
exec root -l "pid/gate_draw_C15d.C(\"pid/${NAME}.json\",\"pid/points_C15d.root\",\"\",\"\",60.0,2.0,${ICLO},${ICHI},$( [[ "$EBEAM" != 0 ]] && echo true || echo false ),${EBEAM},true)"
