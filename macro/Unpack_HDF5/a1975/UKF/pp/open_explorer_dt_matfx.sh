#!/usr/bin/env bash
# Browser explorer for the a1975 16C(d,t)15C MATERIAL/VERTEX test: the SAME dv=1.136 + HDBSCAN
# point clouds fitted twice, differing only in how the triton's energy loss and the vertex
# back-extrapolation are handled:
#
#   matON  : genfit material effects ON  + genfit extrapolation back to the beam axis
#   catima : material effects OFF + hand-applied CATIMA eloss over the vertex gap + back-extrap
#
# Same reco, same PID gate, same geometry, one variable. The page's two-slot fitter switch holds
# the two variants, so the switch and the "overlay the other" checkbox compare them directly.
# Slot labels are rewritten below because make_explorer_html.C always calls the first slot UKF
# and the second GENFIT, and neither word means anything here.
#
#   ./open_explorer_dt_matfx.sh             # Ebeam 180, opens on matON
#   ./open_explorer_dt_matfx.sh 180 catima  # opens on catima instead
#
# The page always activates its FIRST slot on load, so pass the variant you want to land on.
set -o pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EBEAM="${1:-180}"
PRIMARY="${2:-matON}"
SCR="${DTSCR:-/mnt/f/a1975/caches}"   # stable; was a dead session scratchpad
# the cache on disk is dt_kin_maton.root (lower-case "on"); dt_kin_matON.root never existed
CMAT="$SCR/dt_kin_maton.root"; CCAT="$SCR/dt_kin_catima.root"
TMP="${DTTMP:-/tmp}"; OUT="$HOME/a1975_C16_dt_matfx_explorer.html"
MK="$HERE/../../../a2091/UKF/pp/make_explorer_html.C"
for f in "$MK" "$CMAT" "$CCAT"; do [[ -f "$f" ]] || { echo "ERROR: missing $f"; exit 1; }; done

set +u; source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh" >/dev/null 2>&1; set -u

# D2 beam gate [900,1300], not mkexp_pp's H2 default of [950,1350]
root -b -l -q "$HERE/mkexp_pp.C(\"$CMAT\",\"$TMP/exp_dt_matON.root\",1e9,900,1300)"
root -b -l -q "$HERE/mkexp_pp.C(\"$CCAT\",\"$TMP/exp_dt_catima.root\",1e9,900,1300)"
if [[ "$PRIMARY" == "catima" ]]; then
  FIRST="$TMP/exp_dt_catima.root"; SECOND="$TMP/exp_dt_matON.root"; L1="catima"; L2="matFX on"
else
  FIRST="$TMP/exp_dt_matON.root"; SECOND="$TMP/exp_dt_catima.root"; L1="matFX on"; L2="catima"
fi
root -b -l -q "$MK(\"$FIRST\",\"$OUT\",\"16C(d,t)15C  matFX vs CATIMA\",$EBEAM,16.0147013,2.0135532,3.01550072,15.0105993,16,\"0:g.s.,0.740:5/2+,1.218:Sn,3.103,4.220,4.657\",\"$SECOND\")"
[[ -s "$OUT" ]] || { echo "ERROR: explorer not written"; exit 1; }

python3 - "$OUT" "$L1" "$L2" <<'PY'
import sys
p = sys.argv[1]
s = open(p, encoding='utf-8').read()
s = s.replace("const SETS = RAW.ke ? {ukf: RAW} : RAW;",
              "const SETS = RAW.ke ? {ukf: RAW} : RAW;\n"
              "const LBL = {ukf:'%s', genfit:'%s'};\n"
              "const lbl = f => LBL[f] || f;" % (sys.argv[2], sys.argv[3]))
s = s.replace("b.textContent = f;", "b.textContent = lbl(f);")
s = s.replace("${ACTIVE.toUpperCase()}", "${lbl(ACTIVE).toUpperCase()}")
s = s.replace("${other().toUpperCase()} (same cuts)", "${lbl(other()).toUpperCase()} (same cuts)")
s = s.replace("ACTIVE.toUpperCase() + '  |  '", "lbl(ACTIVE).toUpperCase() + '  |  '")
open(p, 'w', encoding='utf-8').write(s)
print("slots relabelled: %s / %s" % (sys.argv[2], sys.argv[3]))
PY

if grep -qi microsoft /proc/version 2>/dev/null; then
  WINHOME="/mnt/c/Users/$(ls /mnt/c/Users | grep -viE 'public|default|all users' | head -1)"
  BROWSER="/mnt/c/Program Files/Google/Chrome/Application/chrome.exe"
  [[ -x "$BROWSER" ]] || BROWSER="/mnt/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"
  cp "$OUT" "$WINHOME/$(basename "$OUT")"
  WINPATH="$(wslpath -w "$WINHOME/$(basename "$OUT")" | sed 's|\\|/|g')"
  echo "opening file:///$WINPATH"
  nohup "$BROWSER" "file:///$WINPATH" >/dev/null 2>&1 &
fi
echo "wrote $OUT"
