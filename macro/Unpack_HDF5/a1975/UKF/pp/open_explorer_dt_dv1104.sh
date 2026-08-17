#!/usr/bin/env bash
# Browser explorer for 16C(d,t)15C at dv 1.10424, Ebeam 184.17 MeV.
#
# NO CUTS ARE PRE-APPLIED. chi2 -> 1e9 and the IC window is opened to 0..1e9, so every
# selection is made in the page, on THIS data. This is a new analysis: nothing is inherited
# from the dv 1.136 / Ebeam 180 work, whose cuts were tuned on a different drift velocity,
# a different track finder and a different z origin.
#
#   ./open_explorer_dt_dv1104.sh [Ebeam]
set -o pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EBEAM="${1:-184.17}"
CACHE="/mnt/f/a1975/caches/dt_kin_dv1104.root"
TMP="${DTTMP:-/mnt/f/a1975/caches/.explorer}"   # stable; the old value was a dead session scratchpad
mkdir -p "$TMP"
OUT="$HOME/a1975_C16_dt_dv1104_explorer.html"
MK="$HERE/../../../a2091/UKF/pp/make_explorer_html.C"
[[ -f "$CACHE" ]] || { echo "ERROR: no cache at $CACHE"; exit 1; }
set +u; source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh" >/dev/null 2>&1; set -u

root -b -l -q "$HERE/mkexp_pp.C(\"$CACHE\",\"$TMP/exp_dt1104.root\",1e9,0,1e9)"
# 15C levels: g.s., 0.740 (5/2+), Sn at 1.218, then the unbound structures.
root -b -l -q "$MK(\"$TMP/exp_dt1104.root\",\"$OUT\",\"16C(d,t)15C  dv 1.10424\",$EBEAM,16.0147013,2.0135532,3.01550072,15.0105993,16,\"0:g.s.,0.740:5/2+,1.218:Sn,3.103,4.220,4.657\",\"\")"
[[ -s "$OUT" ]] || { echo "ERROR: explorer not written"; exit 1; }
python3 - "$OUT" <<'PY'
import sys
p=sys.argv[1]; s=open(p,encoding='utf-8').read()
s=s.replace("const SETS = RAW.ke ? {ukf: RAW} : RAW;",
            "const SETS = RAW.ke ? {ukf: RAW} : RAW;\nconst LBL={ukf:'genfit'};\nconst lbl=f=>LBL[f]||f;")
s=s.replace("b.textContent = f;","b.textContent = lbl(f);")
s=s.replace("${ACTIVE.toUpperCase()}","${lbl(ACTIVE).toUpperCase()}")
open(p,'w',encoding='utf-8').write(s)
PY
python3 "$HERE/add_keoff.py" "$OUT"
cp "$OUT" /mnt/c/Users/Yassid/Desktop/ 2>/dev/null
echo "explorer -> $OUT  (and on the Desktop)"
