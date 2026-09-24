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
# 198 MeV, confirmed. The elastic ridge of pp/ebeam_pp_C15p.C independently returns 199.8-200.4
# over six fit windows (locus estimator 202), so the locus overlay and [Locus check] are drawn on
# a calibrated energy rather than the 195/157/170 guesses that preceded it.
EBEAM="${EBEAM:-198}"
XMAX=45.0
YMAX=1.6
# The plane to draw on. Defaults to the MultiFit+HDBSCAN rebuild of 2026-09-22; set POINTS to
# pid/points_C15p.root to go back to the old PSAMax+TriplClust plane.
POINTS="${POINTS:-pid/points_hdb_C15p.root}"
# Draw the two-body Brho(theta) locus on the plane from the start, not just in [Locus check].
SHOWLOCUS="${SHOWLOCUS:-true}"

set +u
# shellcheck disable=SC1091
source "$REPO/build/config.sh" >/dev/null 2>&1
set -u
[[ -n "${VMCWORKDIR:-}" ]] || { echo "ERROR: config.sh did not set VMCWORKDIR" >&2; exit 1; }
[[ -n "${DISPLAY:-}" ]] || { echo "ERROR: DISPLAY is unset -- the drawer needs an X display" >&2; exit 1; }

cd "$HERE"
[[ -s "$POINTS" ]] || { echo "ERROR: $POINTS missing (run pid/make_points_C15p.C)" >&2; exit 1; }

LOG="/home/yassid/C15p_logs/gate_draw_${species}.log"
mkdir -p "$(dirname "$LOG")"

# Reference gates drawn as overlays (not edited): REFP red, REFD green. Use them to check a
# polygon against the locus of a DIFFERENT channel -- e.g. a gate drawn while the proton locus was
# displayed, re-examined against the (p,d) one.
REFP="${REFP:-}"
REFD="${REFD:-}"
# Output gate file. Defaults to the workspace convention pid/<species>_C15p.json; override to
# keep a different naming (e.g. OUT=pid/deuteron_15Cp.json) without renaming afterwards.
OUT="${OUT:-pid/${species}_C15p.json}"
ARGS="\"${OUT}\",\"${POINTS}\",\"${REFP}\",\"${REFD}\",\
${XMAX},${YMAX},${IC_LO},${IC_HI},${SHOWLOCUS},${EBEAM},true,${Z},${A}"

echo "opening the ${species} gate drawer (Z=${Z} A=${A}) on DISPLAY=${DISPLAY}"
echo "  plane : ${POINTS}, IC [${IC_LO}, ${IC_HI}], single pulse"
echo "  locus : Ebeam ${EBEAM} MeV, overlay ${SHOWLOCUS}  ([Locus check] enabled)"
echo "  bins  : [Locus check] ${LOCUS_NBX:-360} x ${LOCUS_NBY:-240}  (LOCUS_NBX / LOCUS_NBY)"
echo "  out   : ${OUT}   (existing file is backed up to .bak on save)"
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
