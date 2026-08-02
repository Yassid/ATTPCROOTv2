#!/usr/bin/env bash
# Build and open the browser explorer for a1975 16C(d,t)15C (D2 target, triton ejectile).
#
#   ./open_explorer_dt.sh          # Ebeam 180: the resolution study put the g.s. at +0.078 and
#                                  # the 3.103 level at 3.153 with this value, so it stands
#   ./open_explorer_dt.sh 195.5
#
# Cache: the full 47-run ex_dt_a1975.C cache, which unlike the old a1975_panels/dt_kin.root
# carries vertexz and ic -- so the page's z panel works and the D2 beam gate can be applied.
# The gate is [900,1300] (D2 cocktail), NOT the [950,1350] H2 value mkexp_pp defaults to.
#
# Only a GENFIT production exists for this channel, and make_explorer_html.C always puts the first
# cache in the slot it labels "UKF" -- so the label is rewritten here. Without that the page would
# claim genfit tritons were UKF ones.
#
# NOTE: no `set -e`; this ROOT build segfaults in TROOT::EndOfProcessCleanups after writing output.
set -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EBEAM="${1:-180}"
CACHE="${DTCACHE:-/mnt/f/a1975/dt_kin_full.root}"
TMP="${DTTMP:-/tmp}"
mkdir -p "$TMP"
OUT="$HOME/a1975_C16_dt_explorer.html"
MK="$HERE/../../../a2091/UKF/pp/make_explorer_html.C"
[[ -f "$MK" ]] || { echo "ERROR: builder not found at $MK"; exit 1; }
[[ -f "$CACHE" ]] || { echo "ERROR: no (d,t) cache at $CACHE"; exit 1; }

set +u
source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh" >/dev/null 2>&1
set -u

root -b -l -q "$HERE/mkexp_pp.C(\"$CACHE\",\"$TMP/exp_dt.root\",1e9,900,1300)"
# 16C(d,t)15C: beam 16C, target d, ejectile t, residual 15C. 15C levels 0.740/3.103/4.220/4.657.
# The g.s. locus was missing: this list drives the curves drawn on the KE-vs-theta panels, and
# E_x = 0 is the one you actually steer by -- it is the outer edge of the whole kinematic fan,
# so without it there is nothing to check the beam energy or the theta correction against.
root -b -l -q "$MK(\"$TMP/exp_dt.root\",\"$OUT\",\"16C(d,t)15C\",$EBEAM,16.0147013,2.0135532,3.01550072,15.0105993,16,\"0:g.s.,0.740:5/2+,1.218:Sn,3.103,4.220,4.657\",\"\")"
[[ -s "$OUT" ]] || { echo "ERROR: explorer not written"; exit 1; }

python3 - "$OUT" <<'PY'
import sys
p = sys.argv[1]
s = open(p, encoding='utf-8').read()
s = s.replace("const SETS = RAW.ke ? {ukf: RAW} : RAW;",
              "const SETS = RAW.ke ? {ukf: RAW} : RAW;\n"
              "const LBL = {ukf:'genfit'};\n"
              "const lbl = f => LBL[f] || f;")
s = s.replace("b.textContent = f;", "b.textContent = lbl(f);")
s = s.replace("${ACTIVE.toUpperCase()}", "${lbl(ACTIVE).toUpperCase()}")
s = s.replace("${other().toUpperCase()} (same cuts)", "${lbl(other()).toUpperCase()} (same cuts)")
s = s.replace("ACTIVE.toUpperCase() + '  |  '", "lbl(ACTIVE).toUpperCase() + '  |  '")
open(p, 'w', encoding='utf-8').write(s)
print("relabelled: the single slot now reads 'genfit'")
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
