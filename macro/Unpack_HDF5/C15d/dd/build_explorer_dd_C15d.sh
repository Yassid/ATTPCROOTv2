#!/usr/bin/env bash
# Build the (d,d') explorer COMPLETE: bake the cache, then apply every panel and knob.
#
#   dd/build_explorer_dd_C15d.sh [cache.root] [out.html] [Ebeam]
#
# ★ WHY THIS EXISTS. dd/make_explorer_dd_C15d.C produces a BARE page. The IC panel, the per-run Ex
# map, the KE-offset knob and the Ebeam(vz) knob are INJECTED afterwards by four separate scripts,
# and the control defaults by a fifth. Every one is lost on a rebuild, silently, and the page still
# looks fine -- it keeps its numeric icLo/icHi boxes, so nothing appears broken while the IC
# spectrum and both correction knobs are simply absent. dd/open_explorer_dd.sh only ever applied
# two of the five, which is why the (d,d') page has never had the knobs the (p,d) page has.
# Rebuild through here.
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
CACHE="${1:-dd/plots_mfi/dd_kin_C15d.root}"
OUT="${2:-/home/yassid/C15d_dd_explorer_full.html}"
EBEAM="${3:-198.0}"
PP="$HERE/../C15p/pp"          # the injectors live with the C15p page; they take the html as argv
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
cd "$HERE"
[[ -s "$CACHE" ]] || { echo "ERROR: no cache at $CACHE" >&2; exit 1; }

# the IC panel must describe the SAME plane the fits were gated on
case "$CACHE" in
  *plots_mfi*) POINTS=pid/points_mfi_C15d.root; ICJSON=dd/plots_mfi/ic_C15d.json ;;
  *)           POINTS=pid/points_mf_C15d.root;  ICJSON=dd/plots/ic_C15d.json ;;
esac
mkdir -p "$(dirname "$ICJSON")"

echo "=== 1/7  bake the cache into the page ==="
root -b -q "dd/make_explorer_dd_C15d.C(\"$CACHE\",\"$OUT\",\"15C(d,d')\",$EBEAM,15.0105993,2.0141018,2.0141018,15.0105993,15,\"0:g.s.,0.740:5/2+,1.218:Sn\")" \
  2>&1 | grep -E "wrote|tracks|disabled" || true
[[ -s "$OUT" ]] || { echo "ERROR: $OUT not written" >&2; exit 1; }

echo "=== 2/7  IC panel (from $POINTS) ==="
if [[ ! -s "$ICJSON" || "$POINTS" -nt "$ICJSON" ]]; then
  root -b -q "dd/ic_spec_C15d.C(\"$POINTS\",\"$ICJSON\")" 2>&1 | grep -vE '^Info|^Processing|^\s*\||^\s*---|Welcome' | tail -2
fi
[[ -s "$ICJSON" ]] && python3 dd/inject_ic_panel.py "$OUT" "$ICJSON" 2>&1 | tail -1 \
                   || echo "  WARNING: $ICJSON missing -- page will have no IC spectrum"

echo "=== 3/7  per-run Ex map ==="
python3 dd/inject_run_panel.py "$OUT" 2>&1 | tail -1

# ORDER IS LOAD-BEARING: add_keoff introduces the `_k` local that add_ebvz patches against, so
# keoff MUST run first or every ebvz anchor misses and the knob silently does not appear.
echo "=== 4/7  KE-offset knob ==="
python3 "$PP/add_keoff_C15p.py" "$OUT" 2>&1 | tail -2

echo "=== 5/7  Ebeam(vz) knob ==="
python3 "$PP/add_ebvz_C15p.py" "$OUT" 2>&1 | tail -1

# The run panel carries its OWN kine2b, invisible to the two knob injectors, so it keeps plotting
# uncorrected Ex while the main histogram moves. Patch it after both knobs exist.
echo "=== 6/7  per-run map honours the knobs ==="
python3 "$PP/fix_runpanel_corr_C15p.py" "$OUT"

echo "=== 7/7  control defaults ==="
python3 "$PP/set_defaults_C15p.py" "$OUT"

echo; echo "built $OUT  ($(du -h "$OUT" | cut -f1))"
for k in icSpec keoff ebvz; do printf "  %-8s occurrences: %s\n" "$k" "$(grep -c "$k" "$OUT")"; done
