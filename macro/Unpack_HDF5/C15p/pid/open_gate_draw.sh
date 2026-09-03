#!/usr/bin/env bash
# Open the C15p PID gate drawer for one species, detached from the calling shell.
#
#   pid/open_gate_draw.sh deuteron
#   pid/open_gate_draw.sh proton
#   pid/open_gate_draw.sh triton
#
# Two things this handles that a bare `root -l 'pid/gate_draw_C15p.C(...)'` does not:
#
#  1. It calls gApplication->Run() explicitly. The drawer builds a TGMainFrame but never starts
#     an event loop of its own -- interactively that does not matter because the ROOT prompt
#     pumps GUI events, but a detached process has no stdin, ROOT reads EOF and exits, and the
#     window vanishes the instant it appears.
#  2. It bakes in the IC window, so the gate is DRAWN on the same beam-selected plane that
#     apply_gate_C15p.C later applies it to. A gate drawn on the full cocktail and applied to
#     the gated plane selects tracks from a beam it was never meant to include.
#
# The species argument sets Z/A, the output filename and the window title together, so they
# cannot disagree -- they have disagreed before, and a proton gate was saved as a deuteron one.

set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"

species="${1:-deuteron}"
case "$species" in
   proton)   Z=1; A=1 ;;
   deuteron) Z=1; A=2 ;;
   triton)   Z=1; A=3 ;;
   *) echo "unknown species '$species' -- use proton, deuteron or triton" >&2; exit 1 ;;
esac

# IC window from pid/ic_C15p.json (single-pulse) -- the a2091 window, carried over because it was
# chosen on this same detector and beam. Ebeam is what THIS analysis measured
# from its own elastic ridge, used only to draw the optional locus overlay.
IC_LO=979.508
IC_HI=1278.8
EBEAM=202
XMAX=45.0
YMAX=1.6

set +u
# shellcheck disable=SC1091
source "$REPO/build/config.sh" >/dev/null 2>&1
set -u
[[ -n "${VMCWORKDIR:-}" ]] || { echo "ERROR: config.sh did not set VMCWORKDIR" >&2; exit 1; }
[[ -n "${DISPLAY:-}" ]] || { echo "ERROR: DISPLAY is unset -- the drawer needs an X display" >&2; exit 1; }

cd "$HERE"
[[ -s pid/points_C15p.root ]] || { echo "ERROR: pid/points_C15p.root missing (run pid/make_points_C15p.C)" >&2; exit 1; }

LOG="/home/yassid/C15p_logs/gate_draw_${species}.log"
mkdir -p "$(dirname "$LOG")"

ARGS="\"pid/${species}_C15p.json\",\"pid/points_C15p.root\",\"\",\"\",\
${XMAX},${YMAX},${IC_LO},${IC_HI},false,${EBEAM},true,${Z},${A}"

echo "opening the ${species} gate drawer (Z=${Z} A=${A}) on DISPLAY=${DISPLAY}"
echo "  plane : pid/points_C15p.root, IC [${IC_LO}, ${IC_HI}], single pulse"
echo "  out   : pid/${species}_C15p.json   (existing file is backed up to .bak on save)"
echo "  log   : $LOG"

# Invoke as `root -l 'macro.C(args)'`, i.e. cling's .x, and hold stdin open with `tail -f`.
#
# Two things this avoids. Passing the same call through `-e "gROOT->LoadMacro(...); f(...);
# gApplication->Run();"` makes cling compile the whole translation unit eagerly, which dies at
# run time with "symbol '__clang_call_terminate' unresolved while linking" once anything in the
# file touches the C++ stream ABI -- after the plane has been read, so it reads as a GUI crash.
# And redirecting stdin from /dev/null makes ROOT see EOF immediately and exit, taking the
# window with it; the prompt is what pumps the GUI event loop, so it has to stay open.
tail -f /dev/null | setsid nohup root -l "pid/gate_draw_C15p.C(${ARGS})" >"$LOG" 2>&1 &
echo "  pid   : $!"
