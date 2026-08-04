#!/usr/bin/env bash
# Browser explorer for the a1975 16C(d,t)15C DRIFT-VELOCITY test: the same 8 runs reconstructed
# twice, at dv = 1.15 (production) and dv = 1.25, through an identical PID-gated triton chain.
#
# The page's two-slot fitter switch is reused to hold the two drift velocities instead of two
# fitters, so the switch and the "overlay the other" checkbox compare dv directly. The slot
# labels are rewritten below, because make_explorer_html.C always calls the first slot UKF and
# the second GENFIT, and neither word means anything here.
#
#   ./open_explorer_dt_dv.sh            # Ebeam 180, opens on dv 1.15
#   ./open_explorer_dt_dv.sh 191
#   ./open_explorer_dt_dv.sh 180 125    # opens on dv 1.25 instead
#
# The page always activates its FIRST slot on load, so which dv you want to look at has to be
# the one passed first -- hence the second argument, rather than expecting a click.
#
# Measured on these caches: dv 1.25 moves the whole spectrum up 1.20 MeV, taking the g.s. from
# +0.179 to +1.367 and destroying the 3.103 level, so 1.15 is much the better of the two. The
# lever arm, 12.0 MeV per cm/us, puts the g.s. at zero near dv = 1.135.
set -o pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EBEAM="${1:-180}"
PRIMARY="${2:-115}"   # which dv the page opens on
SCR="${DVSCR:-/tmp/claude-1000/-home-yassid/9b2b9ba6-e9b5-4d13-af2f-ad9750f083fc/scratchpad/dv125}"
C115="${SCR}/f8_115g.root"; C125="${SCR}/f8_125.root"
TMP="${DVTMP:-/tmp}"; OUT="$HOME/a1975_C16_dt_dvscan_explorer.html"
MK="$HERE/../../../a2091/UKF/pp/make_explorer_html.C"
for f in "$MK" "$C115" "$C125"; do [[ -f "$f" ]] || { echo "ERROR: missing $f"; exit 1; }; done

set +u; source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh" >/dev/null 2>&1; set -u

# D2 beam gate [900,1300], not mkexp_pp's H2 default of [950,1350]
root -b -l -q "$HERE/mkexp_pp.C(\"$C115\",\"$TMP/exp_dv115.root\",1e9,900,1300)"
root -b -l -q "$HERE/mkexp_pp.C(\"$C125\",\"$TMP/exp_dv125.root\",1e9,900,1300)"
if [[ "$PRIMARY" == "125" ]]; then FIRST="$TMP/exp_dv125.root"; SECOND="$TMP/exp_dv115.root"; L1="dv 1.25"; L2="dv 1.15"; OUT="$HOME/a1975_C16_dt_dv125_explorer.html"
else FIRST="$TMP/exp_dv115.root"; SECOND="$TMP/exp_dv125.root"; L1="dv 1.15"; L2="dv 1.25"; fi
root -b -l -q "$MK(\"$FIRST\",\"$OUT\",\"16C(d,t)15C  dv scan\",$EBEAM,16.0147013,2.0135532,3.01550072,15.0105993,16,\"0:g.s.,0.740:5/2+,1.218:Sn,3.103,4.220,4.657\",\"$SECOND\")"
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
print("slots relabelled: dv 1.15 / dv 1.25")
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
