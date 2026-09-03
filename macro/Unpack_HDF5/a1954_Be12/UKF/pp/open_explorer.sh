#!/usr/bin/env bash
# Rebuild a standalone browser explorer from the kinematics caches and open it in the
# WINDOWS browser. Use this when WSLg/X11 will not show ROOT GUI windows: Chrome runs
# natively on Windows, so nothing here depends on the X server.
#
#   ./open_explorer.sh          # 12Be(p,p')  -- default
#   ./open_explorer.sh pd       # 12Be(p,d)11Be
#   ./open_explorer.sh pt       # 12Be(p,t)10Be
#
# Each page carries BOTH fitters (UKF + GENFIT) and switches between them in-page.
# NOTE: `explorer.exe <file>` silently does nothing on this box -- launch the browser
# binary directly, as below.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CHAN="${1:-pp}"
WINHOME="/mnt/c/Users/$(ls /mnt/c/Users | grep -viE 'public|default|all users' | head -1)"
BROWSER="/mnt/c/Program Files/Google/Chrome/Application/chrome.exe"
[[ -x "$BROWSER" ]] || BROWSER="/mnt/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"

# CACHE is the ex_Be12.C outTag actually on disk: every channel was rebuilt on 2026-08-25 at
# 300 torr + CATIMA with the IC 500-800 beam window, and those caches carry the "800" suffix.
# REF = reference levels of the RESIDUAL nucleus, drawn as kinematic loci on the page.
case "$CHAN" in
  pp) OUT="$HOME/a1954_Be12_pp_explorer.html"; TAG="12Be(p,p')";      MEJ=1.007825; MRES=12.026921
      CACHE=pp800; REF="0:g.s.,2.10:2+_1,4.56:,5.7:" ;;
  pd) OUT="$HOME/a1954_Be12_pd_explorer.html"; TAG="12Be(p,d)11Be";   MEJ=2.014102; MRES=11.021658
      CACHE=pd800; REF="0:g.s. 1/2+,0.320:1/2-,1.778:5/2+,2.654:3/2-,3.400:3/2-" ;;
  pt) OUT="$HOME/a1954_Be12_pt_explorer.html"; TAG="12Be(p,t)10Be";   MEJ=3.016049; MRES=10.013534
      CACHE=pt800; REF="0:g.s. 0+,3.368:2+_1,5.958:2+_2,6.179:0+_2,6.812:Sn,7.371:3-" ;;
  *)  echo "usage: $0 [pp|pd|pt]"; exit 1 ;;
esac

set +u   # thisroot.sh reads unset vars
source "$HOME/fair_install/FairSoft/install/bin/thisroot.sh"
set -u

root -b -l -q "$HERE/make_explorer_html.C(\"$HERE/plots/proton_kin_${CACHE}_ukf.root\",\"$OUT\",\"$TAG\",155,12.026921,1.007825,$MEJ,$MRES,12,\"$REF\",\"$HERE/plots/proton_kin_${CACHE}_genfit.root\")"

# THE ENERGY CORRECTION IS A PATCH, NOT PART OF THE TEMPLATE. add_keoff.py adds the live KE
# offset (applied to the SAME ke in the cut, the theta correction, the maps and the kinematics),
# recolours the reference loci so more than two are visible on a heat map, and gives the
# Ex-vs-theta_cm panel its own zoom. It lives in a1975 and is deliberately NOT copied here --
# one copy, patched for every workspace. A page generated without this step is missing the
# control that says how large an energy bias would have to be to put the states on their lines.
KEOFF="$HERE/../../../a1975/UKF/pp/add_keoff.py"
if [[ -f "$KEOFF" ]]; then
  python3 "$KEOFF" "$OUT"
else
  echo "WARNING: $KEOFF not found -- page written WITHOUT the KE-offset control" >&2
fi

# the browser cannot read \\wsl$ paths reliably -> stage on the Windows side
STAGE="$WINHOME/$(basename "$OUT")"
cp "$OUT" "$STAGE"
[[ -d "$WINHOME/Desktop" ]] && cp "$OUT" "$WINHOME/Desktop/"
WINPATH="$(wslpath -w "$STAGE" | sed 's|\\|/|g')"
echo "opening file:///$WINPATH"
nohup "$BROWSER" "file:///$WINPATH" >/dev/null 2>&1 &
sleep 2
tasklist.exe /FI "IMAGENAME eq $(basename "$BROWSER")" 2>/dev/null | tail -2
