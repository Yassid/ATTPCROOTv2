#!/usr/bin/env bash
# Rebuild the standalone browser explorer from the kinematics caches and open it.
#
#   ./open_explorer.sh          # 15C(p,p')  -- default, gated caches
#   ./open_explorer.sh pd       # 15C(p,d)14C
#   ./open_explorer.sh pp 157   # override the beam energy
#
# The page is a single self-contained HTML file with the data baked in: no server, no ROOT,
# no X11 needed to view it. Each page carries BOTH fitters (UKF + GENFIT no-matFX) and
# switches between them in-page.
#
# This used to be a WSL-only script (it looked for /mnt/c/Users, Windows Chrome, wslpath and
# tasklist.exe) and could not run on this native-Linux box at all. It now uses xdg-open and
# only falls back to the Windows staging dance when actually under WSL.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CHAN="${1:-pp}"
EBEAM="${2:-170}"          # calibrated on the IC+PID-gated sample (was a stale 161 here)

case "$CHAN" in
  pp) OUT="$HOME/a2091_C15_pp_explorer.html"; TAG="15C(p,p')";    MEJ=1.007825; MRES=15.0105993
      # gated caches first, ungated as fallback
      CU="$HERE/plots/proton_kin_g_ukf.root";          [[ -f "$CU" ]] || CU="$HERE/plots/proton_kin_gated.root"
      CG="$HERE/plots/proton_kin_g_genfit_nomat.root"; [[ -f "$CG" ]] || CG=""
      ;;
  # MRES is 14C (14.003242), NOT 13C. The original script passed 13.003355 here, which is 13C:
  # Ex = m4_ex - m4, so the wrong residual offsets the whole spectrum by ~931 MeV.
  pd) OUT="$HOME/a2091_C15_pd_explorer.html"; TAG="15C(p,d)14C";  MEJ=2.014102; MRES=14.003242
      CU="$HERE/plots/proton_kin_pd_ukf.root"
      CG="$HERE/plots/proton_kin_pd_genfit.root";      [[ -f "$CG" ]] || CG=""
      ;;
  *)  echo "usage: $0 [pp|pd] [ebeam]"; exit 1 ;;
esac

if [[ ! -f "$CU" ]]; then
  echo "ERROR: no UKF cache for '$CHAN'. Looked for:"; echo "  $CU"
  echo "Run the Ex step first, e.g.  root -b -q 'pp/ex_C15.C(...)'  to write plots/proton_kin_*.root"
  exit 1
fi
echo "UKF cache    : $CU"
echo "GENFIT cache : ${CG:-<none, single-fitter page>}"

set +u   # thisroot.sh reads unset vars and dies under `set -u`
source "$HOME/fair_install/FairSoft/install/bin/thisroot.sh"
set -u

# beamA = 15 (was 14, left over from the a1954 port)
root -b -l -q "$HERE/make_explorer_html.C(\"$CU\",\"$OUT\",\"$TAG\",$EBEAM,15.0105993,1.007825,$MEJ,$MRES,15,\"\",\"$CG\")"

# ---- open it -------------------------------------------------------------------------
if grep -qi microsoft /proc/version 2>/dev/null; then
  # genuinely under WSL: the Windows browser cannot read \\wsl$ paths reliably, so stage it
  WINHOME="/mnt/c/Users/$(ls /mnt/c/Users | grep -viE 'public|default|all users' | head -1)"
  BROWSER="/mnt/c/Program Files/Google/Chrome/Application/chrome.exe"
  [[ -x "$BROWSER" ]] || BROWSER="/mnt/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"
  cp "$OUT" "$WINHOME/$(basename "$OUT")"
  WINPATH="$(wslpath -w "$WINHOME/$(basename "$OUT")" | sed 's|\\|/|g')"
  echo "opening file:///$WINPATH"
  nohup "$BROWSER" "file:///$WINPATH" >/dev/null 2>&1 &
else
  echo "opening $OUT"
  if command -v xdg-open >/dev/null 2>&1; then nohup xdg-open "$OUT" >/dev/null 2>&1 &
  elif command -v firefox  >/dev/null 2>&1; then nohup firefox  "$OUT" >/dev/null 2>&1 &
  else echo "No browser found. Open this file manually: $OUT"; fi
fi
