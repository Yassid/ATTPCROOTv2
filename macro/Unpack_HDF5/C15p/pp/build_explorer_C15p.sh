#!/usr/bin/env bash
# Build a C15p explorer page COMPLETE: bake the cache, then apply every panel and default.
#
#   pp/build_explorer_C15p.sh [species] [out.html]
#     species: d (deuteron, 15C(p,d)14C -- default) | p (proton, 15C(p,p'))
#
# ★ WHY THIS SCRIPT EXISTS. make_explorer_pp_C15p.C produces a BARE page. The IC panel, the
# per-run Ex map and the Ebeam(vz) knob are INJECTED afterwards by three separate scripts, and the
# sample-hiding control defaults are patched by a fourth. Every one of them is lost on a rebuild,
# silently, and the page still looks fine -- which is how a viewer ends up missing its IC panel and
# plotting a quarter of its tracks with nothing to indicate either. Rebuild through here.
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
SP="${1:-d}"
OUT="${2:-/home/yassid/C15p_hdb_p${SP}_explorer.html}"
CACHE="pp/plots_hdb/pp_kin_C15p_${SP}.root"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
cd "$HERE"
[[ -s "$CACHE" ]] || { echo "ERROR: no cache at $CACHE (run pp/kin_pp_C15p.C with species $SP)" >&2; exit 1; }

if [[ "$SP" == "d" ]]; then
  TAG="15C(p,d)14C"; MEJ=2.01410178; MRES=14.0032420
  REFEX="0:g.s.,6.09:1-,6.59:0+,6.73:3-,7.01:2+,8.18:Sn"      # 14C
else
  TAG="15C(p,p')";   MEJ=1.00782503; MRES=15.0105993
  REFEX="0:g.s.,0.740:5/2+,1.218:Sn"                           # 15C
fi
EBEAM="${EBEAM:-198.0}"

echo "=== 1/5  bake the cache into the page ==="
root -b -q "pp/make_explorer_pp_C15p.C(\"$CACHE\",\"$OUT\",\"$TAG\",$EBEAM,15.0105993,1.00782503,$MEJ,$MRES,15,\"$REFEX\")" \
  2>&1 | grep -E "wrote|tracks|disabled" || true

echo "=== 2/5  IC panel ==="
if [[ -s pid/ic_panel_hdb_C15p.json ]]; then
  python3 pp/inject_ic_panel.py "$OUT" pid/ic_panel_hdb_C15p.json 2>&1 | tail -1
else
  echo "  SKIPPED: pid/ic_panel_hdb_C15p.json missing -- the page will have no IC panel"
fi

echo "=== 3/5  per-run Ex map ==="
python3 pp/inject_run_panel.py "$OUT" 2>&1 | tail -1

echo "=== 4/5  Ebeam(vz) knob ==="
python3 pp/add_ebvz_C15p.py "$OUT" 2>&1 | tail -1

echo "=== 5/5  control defaults (IC filtering OFF -- see set_defaults_C15p.py) ==="
python3 pp/set_defaults_C15p.py "$OUT"

echo
echo "built $OUT  ($(du -h "$OUT" | cut -f1))"
