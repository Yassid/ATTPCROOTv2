#!/usr/bin/env bash
# Rebuild the standalone browser explorer from a kinematics cache and open it in the
# WINDOWS browser. Use this when WSLg/X11 will not show ROOT GUI windows: Chrome runs
# natively on Windows, so nothing here depends on the X server.
#
#   ./open_explorer.sh                              # clean155 cache, 155 MeV
#   ./open_explorer.sh plots/proton_kin_pd.root      # another cache
#
# NOTE: `explorer.exe <file>` silently does nothing on this box -- launch the browser
# binary directly, as below.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CACHE="${1:-}"
OUT="$HOME/a1954_Be12_explorer.html"
WINDIR="/mnt/c/Users/$(basename "$(ls -d /mnt/c/Users/* | grep -viE 'public|default|all users' | head -1)")"
BROWSER="/mnt/c/Program Files/Google/Chrome/Application/chrome.exe"
[[ -x "$BROWSER" ]] || BROWSER="/mnt/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"

set +u   # thisroot.sh reads unset vars
source "$HOME/fair_install/FairSoft/install/bin/thisroot.sh"
set -u
if [[ -n "$CACHE" ]]; then
   root -b -l -q "$HERE/make_explorer_html.C(\"$CACHE\",\"$OUT\")"
else
   root -b -l -q "$HERE/make_explorer_html.C()"
fi

# the browser cannot read \\wsl$ paths reliably -> stage on the Windows side
STAGE="$WINDIR/a1954_Be12_explorer.html"
cp "$OUT" "$STAGE"
WINPATH="$(wslpath -w "$STAGE" | sed 's|\\|/|g')"
echo "opening file:///$WINPATH"
nohup "$BROWSER" "file:///$WINPATH" >/dev/null 2>&1 &
sleep 2
tasklist.exe /FI "IMAGENAME eq $(basename "$BROWSER")" 2>/dev/null | tail -2
