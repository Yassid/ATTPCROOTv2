#!/usr/bin/env bash
# Build and open the browser explorer for a1975 16C(p,d)15C.
#
#   ./open_explorer_pd.sh              # tight vs loose gate, Ebeam 195.5
#   ./open_explorer_pd.sh 192          # override the beam energy
#
# The two in-page "fitter" slots are used here for the two GATE versions -- same fitter (genfit),
# same runs, same everything else, so the switch isolates the deuteron PID gate. a2091's
# make_explorer_html.C hardcodes the slot names 'ukf'/'genfit' and renders the buttons straight
# from those keys, so the labels are rewritten afterwards; without that the page would claim to be
# comparing fitters when it is comparing gates.
#
# NOTE: no `set -e` -- this ROOT build segfaults in TROOT::EndOfProcessCleanups AFTER writing its
# output, so every root call returns non-zero on success.
set -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CACHE="${PDCACHE:-/tmp/claude-1000/-home-yassid/0335450a-a4bc-427a-9507-0b3184e093bf/scratchpad/pdcache}"
EBEAM="${1:-195.5}"
OUT="$HOME/a1975_C16_pd_explorer.html"
MK="$HERE/../../../a2091/UKF/pp/make_explorer_html.C"
[[ -f "$MK" ]] || { echo "ERROR: explorer builder not found at $MK"; exit 1; }

set +u
source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh" >/dev/null 2>&1
set -u

T_IN="$CACHE/pd_kin_tight.root"; L_IN="$CACHE/pd_kin_loose.root"
[[ -f "$T_IN" ]] || { echo "ERROR: no tight cache at $T_IN"; exit 1; }
echo "tight gate : $T_IN"
root -b -l -q "$HERE/mkexp_pp.C(\"$T_IN\",\"$CACHE/exp_tight.root\")"
CU="$CACHE/exp_tight.root"; CG=""
if [[ -f "$L_IN" ]]; then
  echo "loose gate : $L_IN"
  root -b -l -q "$HERE/mkexp_pp.C(\"$L_IN\",\"$CACHE/exp_loose.root\")"
  CG="$CACHE/exp_loose.root"
else
  echo "WARNING: no loose cache -- single-gate page."
fi

# 16C(p,d)15C: beam 16C, target p, ejectile d, residual 15C, beamA=16.
# 15C reference levels: g.s., 0.740 (5/2+), 3.103, 4.220, 4.657 MeV.
root -b -l -q "$MK(\"$CU\",\"$OUT\",\"16C(p,d)15C\",$EBEAM,16.0147013,1.00782503,2.0135532,15.0105993,16,\"0.740,3.103,4.220,4.657\",\"$CG\")"
[[ -s "$OUT" ]] || { echo "ERROR: explorer not written"; exit 1; }

# ---- relabel the two slots: they hold GATES here, not fitters ----
python3 - "$OUT" <<'PY'
import sys, re
p = sys.argv[1]
s = open(p, encoding='utf-8').read()
# a label map consulted wherever the raw slot key would otherwise be shown
s = s.replace("const SETS = RAW.ke ? {ukf: RAW} : RAW;",
              "const SETS = RAW.ke ? {ukf: RAW} : RAW;\n"
              "const LBL = {ukf:'tight gate', genfit:'loose gate'};\n"
              "const lbl = f => LBL[f] || f;")
s = s.replace("b.textContent = f;", "b.textContent = lbl(f);")
s = s.replace("${ACTIVE.toUpperCase()}", "${lbl(ACTIVE).toUpperCase()}")
s = s.replace("${other().toUpperCase()} (same cuts)", "${lbl(other()).toUpperCase()} (same cuts)")
s = s.replace("ACTIVE.toUpperCase() + '  |  '", "lbl(ACTIVE).toUpperCase() + '  |  '")
s = s.replace("overlay the other fitter", "overlay the other gate")
s = s.replace(">fitter<", ">deuteron PID gate<")
open(p, 'w', encoding='utf-8').write(s)
print("relabelled: slots now read 'tight gate' / 'loose gate'")
PY

if grep -qi microsoft /proc/version 2>/dev/null; then
  WINHOME="/mnt/c/Users/$(ls /mnt/c/Users | grep -viE 'public|default|all users' | head -1)"
  BROWSER="/mnt/c/Program Files/Google/Chrome/Application/chrome.exe"
  [[ -x "$BROWSER" ]] || BROWSER="/mnt/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"
  cp "$OUT" "$WINHOME/$(basename "$OUT")"
  WINPATH="$(wslpath -w "$WINHOME/$(basename "$OUT")" | sed 's|\\|/|g')"
  echo "opening file:///$WINPATH"
  nohup "$BROWSER" "file:///$WINPATH" >/dev/null 2>&1 &
else
  command -v xdg-open >/dev/null 2>&1 && nohup xdg-open "$OUT" >/dev/null 2>&1 || echo "open manually: $OUT"
fi
echo "wrote $OUT"
