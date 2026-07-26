#!/usr/bin/env bash
# Rebuild a standalone browser explorer from the kinematics caches and open it in the
# WINDOWS browser. Use this when WSLg/X11 will not show ROOT GUI windows: Chrome runs
# natively on Windows, so nothing here depends on the X server.
#
#   ./open_explorer.sh          # 14C(p,p')  -- default
#   ./open_explorer.sh pd       # 14C(p,d)13C
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

case "$CHAN" in
  pp) OUT="$HOME/a1954_C14_pp_explorer.html"; TAG="14C(p,p')";      MEJ=1.007825; MRES=14.003242 ;;
  pd) OUT="$HOME/a1954_C14_pd_explorer.html"; TAG="14C(p,d)13C";   MEJ=2.014102; MRES=13.003355 ;;
  *)  echo "usage: $0 [pp|pd]"; exit 1 ;;
esac

set +u   # thisroot.sh reads unset vars
source "$HOME/fair_install/FairSoft/install/bin/thisroot.sh"
set -u

root -b -l -q "$HERE/make_explorer_html.C(\"$HERE/plots/proton_kin_${CHAN}_ukf.root\",\"$OUT\",\"$TAG\",161,14.003242,1.007825,$MEJ,$MRES,14,\"\",\"$HERE/plots/proton_kin_${CHAN}_genfit.root\")"

# the browser cannot read \\wsl$ paths reliably -> stage on the Windows side
STAGE="$WINHOME/$(basename "$OUT")"
cp "$OUT" "$STAGE"
[[ -d "$WINHOME/Desktop" ]] && cp "$OUT" "$WINHOME/Desktop/"
WINPATH="$(wslpath -w "$STAGE" | sed 's|\\|/|g')"
echo "opening file:///$WINPATH"
nohup "$BROWSER" "file:///$WINPATH" >/dev/null 2>&1 &
sleep 2
tasklist.exe /FI "IMAGENAME eq $(basename "$BROWSER")" 2>/dev/null | tail -2
