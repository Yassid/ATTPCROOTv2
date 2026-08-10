#!/usr/bin/env bash
# Build and open the standalone browser explorer for a1975 16C(p,p').
#
#   ./open_explorer_pp.sh                 # both fitters, Ebeam 195.5 (spacing-calibrated)
#   ./open_explorer_pp.sh 192             # override the beam energy
#
# Reuses a2091's make_explorer_html.C -- the page is self-contained (data baked in), so no
# server, no ROOT and no X11 are needed to view it. Ebeam, cuts, binning and the theta
# correction are all live in-page.
# NOTE: no `set -e`. This ROOT build segfaults in TROOT::EndOfProcessCleanups AFTER the macro
# has finished and written its output, so every root call returns non-zero even on success;
# with `set -e` the script would abort after the first one. Outputs are checked explicitly.
set -o pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CACHE="${PPCACHE:-/tmp/claude-1000/-home-yassid/0335450a-a4bc-427a-9507-0b3184e093bf/scratchpad/ppcache}"
EBEAM="${1:-195.5}"
OUT="$HOME/a1975_C16_pp_explorer.html"
# $HERE is .../a1975/UKF/pp, so a2091 is THREE levels up (not two -- that pointed inside a1975)
MK="$HERE/../../../a2091/UKF/pp/make_explorer_html.C"
[[ -f "$MK" ]] || { echo "ERROR: explorer builder not found at $MK"; exit 1; }

set +u
source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh" >/dev/null 2>&1
set -u

# genfit is the reference cache; UKF is optional (its per-run pass may still be building)
GF_IN="$CACHE/pp_kin_genfit.root"; [[ -f "$GF_IN" ]] || GF_IN="$CACHE/partial.root"
# prefer the PID-gated UKF cache (pp/cache_pp_gated.C): the genfit production was gated in-fitter,
# so the ungated pp_kin_ukf.root would compare ~1.9x more tracks against it and is not like-for-like
UK_IN="$CACHE/pp_kin_ukf_gated.root"
if [[ ! -f "$UK_IN" ]]; then
  UK_IN="$CACHE/pp_kin_ukf.root"
  [[ -f "$UK_IN" ]] && echo "WARNING: using the UNGATED UKF cache -- selections do not match genfit."
fi
[[ -f "$GF_IN" ]] || { echo "ERROR: no genfit cache at $GF_IN"; exit 1; }

echo "genfit cache : $GF_IN"
root -b -l -q "$HERE/mkexp_pp.C(\"$GF_IN\",\"$CACHE/exp_genfit.root\")"
CG="$CACHE/exp_genfit.root"
if [[ -f "$UK_IN" ]]; then
  echo "UKF cache    : $UK_IN"
  root -b -l -q "$HERE/mkexp_pp.C(\"$UK_IN\",\"$CACHE/exp_ukf.root\")"
  CU="$CACHE/exp_ukf.root"
else
  # make_explorer_html.C hardcodes the two button labels, so a single cache always lands in
  # the slot labelled "UKF". Say so plainly rather than letting the page mislabel silently.
  echo "UKF cache    : <not built yet>"
  echo "WARNING: single-fitter page -- the button labelled UKF is showing GENFIT data."
  CU="$CG"; CG=""
fi

# 16C(p,p'): beam 16C, target p, ejectile p, residual 16C, beamA=16.
# Reference levels. THE GROUND STATE MUST BE IN THIS LIST: the builder draws one kinematic curve
# per entry and nothing else, so omitting 0 leaves the elastic locus -- the strongest thing in the
# data and the one the beam energy is tuned against -- with no line to compare to. It was missing
# here, which made the page look as though the g.s. curve had failed to draw.
# 16C: 0+ g.s., 2+ 1.766, 0+ 3.027, 2+ 3.986, 4+ 4.142 MeV.
root -b -l -q "$MK(\"$CU\",\"$OUT\",\"16C(p,p')\",$EBEAM,16.0147013,1.00782503,1.00782503,16.0147013,16,\"0:g.s.,1.766:2+,3.027:0+,3.986:2+,4.142:4+\",\"$CG\")"

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
