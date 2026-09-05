#!/usr/bin/env bash
# Build and open the 46Ar(3He,d)47K kinematics explorer.
#
#   ./explorer/open_explorer_Ar46.sh              # rebuild the page and open it
#   NOOPEN=1 ./explorer/open_explorer_Ar46.sh     # rebuild only
#
# THREE STEPS, AND THE MIDDLE ONE IS NOT OPTIONAL:
#   kin_Ar46.C            per configuration -> plots/kin_<tag>.root   (skipped if present)
#   make_explorer_Ar46.C  N caches + the SHARED template -> one page
#   patch_ar46.py         the 46Ar beam: per-event E_beam(z) and a 598 MeV slider
#
# Without the patch the page still loads and still draws -- with every kinematic locus at one
# constant beam energy the data never actually has, since 46Ar loses 95.7 MeV crossing the drift.
# That is a page that lies rather than one that fails, so this script refuses to open an
# unpatched file.
#
# THE BROWSER CANNOT READ \\wsl$ PATHS RELIABLY, so the page is staged on the Windows side first.
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DIR="$(dirname "$HERE")"
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
OUT="${OUT:-$HOME/ar46_3Hed_explorer.html}"

# tag : fit directory : sim directory (generation is shared by the two pad arms of a field)
A=/mnt/f/ar46_3hed_OLD_2.85T_placeholder
B=/mnt/f/ar46_3hed_gen_B39
CFG=(
  "b285_attpc:/mnt/f/ar46_3hed_mx_B285_attpc:$A"
  "b39_attpc:/mnt/f/ar46_3hed_mx_B39_attpc:$B"
  "b285_2mm:/mnt/f/ar46_3hed_mx_B285_2mm:$A"
  "b39_2mm:/mnt/f/ar46_3hed_mx_B39_2mm:$B"
)

set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export VMCWORKDIR="$REPO"
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
cd "$DIR"

for c in "${CFG[@]}"; do
  IFS=: read -r tag fit sim <<< "$c"
  # REBUILD ONLY WHAT IS MISSING, and test the PRODUCT rather than its name: a cache that exists
  # but holds no pk tree is exactly what a killed job leaves behind.
  n=$(root -b -q -l -e "TFile*f=TFile::Open(\"plots/kin_${tag}.root\");TTree*t=(f&&!f->IsZombie())?(TTree*)f->Get(\"pk\"):nullptr;printf(\"N %lld\\n\",t?t->GetEntries():-1);" 2>/dev/null | awk '/^N /{print $2}' || true)
  if [ "${n:--1}" -gt 0 ]; then
    echo "[cache] $tag already built ($n tracks)"
  else
    echo "[cache] $tag building from $fit"
    root -b -q -l "kin_Ar46.C(\"$fit\",\"$sim\",\"$tag\")" || { echo "kin_Ar46 FAILED for $tag"; exit 1; }
  fi
done

# Build the set list from CFG so the launcher and the page cannot disagree about which
# configurations exist -- and so adding an arm here is the only edit needed.
SETS=""
for c in "${CFG[@]}"; do
  IFS=: read -r tag fit sim <<< "$c"
  case "$tag" in
    b285_attpc) nm="2.85T AT-TPC" ;; b39_attpc) nm="3.9T AT-TPC" ;;
    b285_2mm)   nm="2.85T 2mm"    ;; b39_2mm)   nm="3.9T 2mm"    ;;
    *)          nm="$tag"         ;;
  esac
  SETS="${SETS:+$SETS,}${nm}=../plots/kin_${tag}.root"
done
echo "[sets] $SETS"

root -b -q -l "explorer/make_explorer_Ar46.C(\"$SETS\",\"$OUT\")" || { echo "generator FAILED"; exit 1; }
python3 "$HERE/patch_ar46.py" "$OUT"
grep -q "__AR46_PATCHED__" "$OUT" || { echo "REFUSING TO OPEN: $OUT is not patched"; exit 1; }

[ -n "${NOOPEN:-}" ] && { echo "built $OUT (not opened)"; exit 0; }

WINHOME="/mnt/c/Users/$(ls /mnt/c/Users | grep -viE 'public|default|all users' | head -1)"
# STAGE ON THE DESKTOP so the page is where collaborators are told to look, and pick the LIVE
# Desktop: OneDrive redirection leaves two of them and the stale one still exists. Choose by
# most-recent content, not by which path exists -- both do.
DESK="$WINHOME/Desktop"
if [ -d "$WINHOME/OneDrive/Desktop" ]; then
  a=$(find "$WINHOME/Desktop" -maxdepth 1 -printf '%T@\n' 2>/dev/null | sort -n | tail -1)
  b=$(find "$WINHOME/OneDrive/Desktop" -maxdepth 1 -printf '%T@\n' 2>/dev/null | sort -n | tail -1)
  awk -v x="${a:-0}" -v y="${b:-0}" 'BEGIN{exit !(y>x)}' && DESK="$WINHOME/OneDrive/Desktop"
fi
mkdir -p "$DESK"
STAGE="$DESK/$(basename "$OUT")"
cp "$OUT" "$STAGE"
echo "[desktop] $STAGE"
WINPATH=$(printf '%s' "$STAGE" | sed 's|/mnt/c|C:|; s|/|\\|g')
BROWSER="/mnt/c/Program Files/Google/Chrome/Application/chrome.exe"
[[ -x "$BROWSER" ]] || BROWSER="/mnt/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"
if [[ -x "$BROWSER" ]]; then
  "$BROWSER" "$WINPATH" >/dev/null 2>&1 &
  echo "opened $WINPATH"
else
  # `explorer.exe <file>` silently does nothing on this box; cmd start is the reliable route.
  cmd.exe /c start "" "$WINPATH" >/dev/null 2>&1 || true
  echo "staged $WINPATH -- open it by hand if the browser did not appear"
fi
