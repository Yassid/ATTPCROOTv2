#!/usr/bin/env bash
# Open the ion-chamber (beam) gate picker.
#
#   ./pid/open_icgate_gui.sh                    # pick the window with two clicks
#   ./pid/open_icgate_gui.sh ic_C15d 900 1792   # or set it directly
#
# LEFT-CLICK once at the LOW edge and once at the HIGH edge; only the x positions are used.
# Writes pid/<name>.json.
#
# Interactive: needs a display, and ROOT's prompt quits on stdin EOF, so stdin is held open when
# this is launched detached.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
NAME="${1:-ic_C15d}"; LO="${2:--1}"; HI="${3:--1}"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
command -v root >/dev/null || { echo "ERROR: root not on PATH" >&2; exit 1; }
[[ -n "${DISPLAY:-}${WAYLAND_DISPLAY:-}" ]] || { echo "ERROR: no display" >&2; exit 1; }
cd "$HERE"; mkdir -p pid/plots
[[ -t 0 ]] || exec 0< <(while :; do sleep 3600; done)
echo "IC gate picker: name=$NAME window=[$LO,$HI]  (two clicks if unset)"
exec root -l "pid/draw_icgate_C15d.C(\"$NAME\",$LO,$HI)"
