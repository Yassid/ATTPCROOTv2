#!/usr/bin/env bash
# Open the raw-vs-matched PID comparison GUI, detached.
#
#   pid/open_pid_compare.sh
#
# Two things this handles, both learned the hard way with the gate drawer:
#  * `root -l 'macro.C(args)'` (cling's .x), NOT `-e "gROOT->LoadMacro(...)"`. The -e form compiles
#    the whole translation unit eagerly and dies at run time with "symbol '__clang_call_terminate'
#    unresolved" once anything touches the C++ stream ABI -- after the data is read, so it reads
#    as a GUI crash rather than a load failure.
#  * stdin held open with `tail -f`. Redirecting from /dev/null makes ROOT see EOF and exit
#    immediately, taking the window with it; the prompt is what pumps the GUI event loop.
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
[[ -n "${DISPLAY:-}" ]] || { echo "ERROR: DISPLAY is unset -- this needs an X display" >&2; exit 1; }
cd "$HERE"
[[ -n "$(ls -1 /home/yassid/C15d_reco/*_pid.root 2>/dev/null)" ]] \
  || { echo "ERROR: no _pid.root ntuples in /home/yassid/C15d_reco" >&2; exit 1; }
LOG=/home/yassid/C15d_logs/pid_compare.log
mkdir -p "$(dirname "$LOG")"
echo "opening the PID comparison GUI on DISPLAY=$DISPLAY"
echo "  runs : $(ls -1 /home/yassid/C15d_reco/*_pid.root | wc -l) pid ntuples"
echo "  gain : gainmatch_C15d.csv"
echo "  log  : $LOG"
tail -f /dev/null | setsid nohup root -l "pid_compare_C15d.C()" >"$LOG" 2>&1 &
echo "  pid  : $!"
