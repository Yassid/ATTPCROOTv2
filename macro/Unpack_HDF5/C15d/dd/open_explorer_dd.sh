#!/usr/bin/env bash
# Rebuild the (d,d') kinematics explorer from the gated deuteron ntuple and open it.
#
#   ./dd/open_explorer_dd.sh            # Ebeam 90 MeV, from this analysis's own elastic ridge
#   ./dd/open_explorer_dd.sh 95         # try another beam energy
#
# One self-contained HTML file: Ex and theta_cm are recomputed in the page from (KE, theta_lab)
# with the same two-body expressions, so the beam energy and every cut stay live in the browser.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
EBEAM="${1:-90}"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
cd "$HERE"
[[ -s dd/plots/dd_kin_C15d.root ]] || { echo "ERROR: run dd/kin_dd_C15d.C first" >&2; exit 1; }
OUT="$HOME/C15d_dd_explorer.html"
root -b -q "dd/make_explorer_dd_C15d.C(\"\",\"$OUT\",\"15C(d,d')\",$EBEAM)" 2>&1 | grep -vE '^Info|^Processing|^\s*\||^\s*---|Welcome'
[[ -s "$OUT" ]] || { echo "ERROR: $OUT not written" >&2; exit 1; }
if [[ -n "${DISPLAY:-}${WAYLAND_DISPLAY:-}" ]]; then
   setsid xdg-open "$OUT" >/dev/null 2>&1 </dev/null & disown || true
   echo "opened $OUT"
else
   echo "no display -- open $OUT yourself"
fi
