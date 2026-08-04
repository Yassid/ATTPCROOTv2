#!/usr/bin/env bash
# Browser explorer for the a1975 16C(d,t)15C DRIFT-VELOCITY test: two full 47-run productions
# that differ ONLY in DriftVelocity. Everything else is held fixed -- same hdbscan reco config
# (thr 40, mcs 0/3), same triton PID gate, same geometry, and the same fit configuration
# (genfit matEffects OFF + back-extrapolation + hand-applied CATIMA over the vertex gap).
#
#   dv 1.136 : reco_d2_hdb1136 -> gf_dt_catima
#   dv 1.35  : reco_d2_dv135   -> gf_dt_dv135
#
# Pairing 1.35 against gf_dt_matON instead would compare drift velocity AND material treatment
# at once, which is how an earlier comparison came out backwards. Use the catima cache.
#
# Ebeam defaults to 184 MeV = 11.5 MeV/u, the correct beam energy.
# The page always activates its FIRST slot on load, so pass the one you want to land on.
#
#   ./open_explorer_dt_dvcmp.sh            # Ebeam 184, opens on dv 1.136
#   ./open_explorer_dt_dvcmp.sh 184 135    # opens on dv 1.35 instead
#   ./open_explorer_dt_dvcmp.sh 184 1136 135    # explicit pair
#   ./open_explorer_dt_dvcmp.sh 184 1136 150    # 1.136 against 1.50
# The FIRST tag is the slot the page opens on.
set -o pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EBEAM="${1:-184}"
TAGA="${2:-1136}"; TAGB="${3:-135}"
CDIR="${DTCACHE:-/mnt/f/a1975/caches}"
cachefor(){ case "$1" in
    1136) echo "$CDIR/dt_kin_catima_dv1136.root";;
    135)  echo "$CDIR/dt_kin_dv135.root";;
    150)  echo "$CDIR/dt_kin_dv150.root";;
    *)    echo "";; esac; }
lblfor(){ case "$1" in 1136) echo "dv 1.136";; 135) echo "dv 1.35";; 150) echo "dv 1.50";; *) echo "dv $1";; esac; }
CA="$(cachefor "$TAGA")"; CB="$(cachefor "$TAGB")"
[[ -n "$CA" && -n "$CB" ]] || { echo "ERROR: unknown dv tag (use 1136 | 135 | 150)"; exit 1; }
# Output is named after the PAIR: the page holds two slots (its switch toggles between exactly
# two), so three drift velocities means three pages rather than one three-way page.
TMP="${DTTMP:-/tmp}"; OUT="$HOME/a1975_C16_dt_dv${TAGA}_vs_${TAGB}_explorer.html"
MK="$HERE/../../../a2091/UKF/pp/make_explorer_html.C"
for f in "$MK" "$CA" "$CB"; do [[ -f "$f" ]] || { echo "ERROR: missing $f"; exit 1; }; done

set +u; source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$HOME/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}"

# D2 beam gate [900,1300], not mkexp_pp's H2 default of [950,1350]
root -b -l -q "$HERE/mkexp_pp.C(\"$CA\",\"$TMP/exp_dt_${TAGA}.root\",1e9,900,1300)"
root -b -l -q "$HERE/mkexp_pp.C(\"$CB\",\"$TMP/exp_dt_${TAGB}.root\",1e9,900,1300)"
FIRST="$TMP/exp_dt_${TAGA}.root"; SECOND="$TMP/exp_dt_${TAGB}.root"
L1="$(lblfor "$TAGA")"; L2="$(lblfor "$TAGB")"
root -b -l -q "$MK(\"$FIRST\",\"$OUT\",\"16C(d,t)15C  $L1 vs $L2\",$EBEAM,16.0147013,2.0135532,3.01550072,15.0105993,16,\"0:g.s.,0.740:5/2+,1.218:Sn,3.103,4.220,4.657\",\"$SECOND\")"
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
  cp "$OUT" "$WINHOME/Desktop/$(basename "$OUT")" 2>/dev/null || cp "$OUT" "$WINHOME/$(basename "$OUT")"
  WINPATH="$(wslpath -w "$WINHOME/Desktop/$(basename "$OUT")" 2>/dev/null | sed 's|\\|/|g')"
  echo "opening file:///$WINPATH"
  nohup "$BROWSER" "file:///$WINPATH" >/dev/null 2>&1 &
fi
echo "wrote $OUT"
