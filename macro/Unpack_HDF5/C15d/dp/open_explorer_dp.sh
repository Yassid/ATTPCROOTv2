#!/usr/bin/env bash
# Rebuild the (d,p) kinematics explorer from the gated proton ntuple and open it.
#
#   ./dp/open_explorer_dp.sh            # Ebeam 207 MeV, from the (d,d) elastic ridge
#   ./dp/open_explorer_dp.sh 205        # try another beam energy
#
# ★ 15C(d,p)16C: the ejectile is the PROTON and the residual is 16C. Unlike (d,d'), this channel
# is NOT degenerate in beam mass -- 16C's level scheme (2+ at 1.766) is nothing like 15C's.
#
# One self-contained HTML file: Ex and theta_cm are recomputed in the page from (KE, theta_lab)
# with the same two-body expressions, so the beam energy and every cut stay live in the browser.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
EBEAM="${1:-207}"
# 16C levels, so the g.s. band can be told from its neighbours. These are the RESIDUAL's levels --
# the whole reason (d,p) discriminates the beam where (d,d') cannot.
REFEX="0:g.s.,1.766:2+,3.027:0+,3.986:2+,4.25:Sn"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
cd "$HERE"
[[ -s dp/plots/dp_kin_C15d.root ]] || { echo "ERROR: run dp/kin_dp_C15d.C first" >&2; exit 1; }
OUT="$HOME/C15d_dp_explorer.html"
root -b -q "dd/make_explorer_dd_C15d.C(\"dp/plots/dp_kin_C15d.root\",\"$OUT\",\"15C(d,p)\",$EBEAM,15.0105993,2.0141018,1.00782503207,16.0147413,15,\"$REFEX\")" 2>&1 | grep -vE '^Info|^Processing|^\s*\||^\s*---|Welcome'
[[ -s "$OUT" ]] || { echo "ERROR: $OUT not written" >&2; exit 1; }

# ★ The ion-chamber panel is INJECTED, not part of the shared template, so it has to be re-applied
# on every regeneration. It was not wired in here before and was silently lost for six days: the
# page kept its numeric icLo/icHi boxes, so nothing looked broken -- only the spectrum was gone.
# The dump is cheap (one pass over points_C15d.root) but not free, so it is rebuilt only when the
# points file is newer than the JSON.
IC_JSON="dp/plots/ic_C15d.json"
if [[ ! -s "$IC_JSON" || pid/points_C15d.root -nt "$IC_JSON" ]]; then
   root -b -q 'dd/ic_spec_C15d.C("pid/points_C15d.root","dp/plots/ic_C15d.json")' 2>&1 | grep -vE '^Info|^Processing|^\s*\||^\s*---|Welcome'
fi
if [[ -s "$IC_JSON" ]]; then
   python3 dd/inject_ic_panel.py "$OUT" "$IC_JSON"
else
   echo "WARNING: $IC_JSON missing -- page has the IC boxes but no spectrum panel" >&2
fi

# The per-run Ex map is injected the same way, and was lost the same way. It needs no dump of its
# own -- it calls the page's own kine2b/thCorr/pass on the run column already in the cache.
python3 dd/inject_run_panel.py "$OUT"
if [[ -n "${DISPLAY:-}${WAYLAND_DISPLAY:-}" ]]; then
   setsid xdg-open "$OUT" >/dev/null 2>&1 </dev/null & disown || true
   echo "opened $OUT"
else
   echo "no display -- open $OUT yourself"
fi
