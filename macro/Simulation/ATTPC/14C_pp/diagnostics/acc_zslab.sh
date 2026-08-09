#!/usr/bin/env bash
# Recompute the GENFIT acceptance inside a vertex-z slab, for every level, and merge the seeds.
#
# WHY A SLAB. The excited states only populate the first ~40 cm of the chamber while the elastic
# fills all 95 cm, so a yield ratio between them is a ratio over two different target thicknesses
# unless both are restricted to the same z window. The absolute normalisation is built by scaling
# the excited-state yields to the elastic one, so this matters directly: get it wrong and every
# deformation length is wrong by the thickness ratio.
#
# The window also happens to be where the vertex is measurable at all -- reco z tracks truth to
# about 15 mm below 600 mm, but at a truth z of 890 mm the mean reconstructed z is 409 mm.
#
# No re-simulation: this re-reads the existing sim + genfit files, so it is minutes not hours.
#
# chi2Cut defaults to 1e9 (no cut), and that is almost certainly the one to use. The simulation
# does not reproduce the fit-quality distribution of the data at all: chi2/ndf < 5 keeps ~93 % of
# simulated events but only 65 % of real ones, and in data the survival runs from 45 % at
# theta_cm 20-40 to 92 % at 100-120. A simulated acceptance cannot correct a cut it does not
# reproduce, so the cut acceptance is provided only for consistency studies.
#
#   ./acc_zslab.sh [zMin] [zMax] [chi2Cut] [outSuffix]
#
# chi2Cut MUST MATCH THE DATA CACHE. proton_kin_300gfx_ex.root has chi2/ndf < 5 applied while
# proton_kin_300gfx_nc.root has no cut at all; correcting the first with a no-chi2 acceptance
# credits the data with tracks the selection threw away, and for an ABSOLUTE normalisation that
# is a direct multiplicative error rather than a shape nuisance.
set -eo pipefail
ZMIN=${1:-10}; ZMAX=${2:-400}; CHI2=${3:-1e9}; SUF=${4:-}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); SIM=$(cd "$HERE/.." && pwd)
REPO=$(cd "$SIM/../../../.." && pwd)
ACC=/mnt/f/a1954_C14_acc_gf
OUT=/mnt/f/a1954_C14_acc_gf_z${ZMIN}_${ZMAX}${SUF}
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
cd "$SIM"; mkdir -p "$OUT"

run_level () {
  local TAG=$1 EX=$2 SEEDS=$3
  local used=""
  for s in $SEEDS; do
    local J="${TAG}_s${s}"
    [ -f "$ACC/$J.marker" ] || { echo "  $J: no marker, skipping"; continue; }
    root -b -q -l "acceptance_C14.C(\"$ACC/${J}_sim.root\",\"$ACC/${J}_genfit.root\",\"$J\",$EX,161.0,$CHI2,36,180.0,10.0,0.5,2.0,kTRUE,$ZMIN,$ZMAX)" \
         > "$OUT/${J}_acc.log" 2>&1
    if grep -q "overall acceptance" "$OUT/${J}_acc.log"; then
      mv -f "acceptance_${J}.root" "$OUT/" 2>/dev/null || mv -f "$HERE/acceptance_${J}.root" "$OUT/" 2>/dev/null || true
      rm -f "acceptance_${J}.png" "$HERE/acceptance_${J}.png"
      used="${used}${used:+,}${s}"
      echo "  $J: $(grep 'overall acceptance' "$OUT/${J}_acc.log")"
    else
      echo "  $J: ACC_FAILED -- see $OUT/${J}_acc.log"
    fi
  done
  # Merging "whatever worked" without saying so would hide a missing seed in the error bar.
  [ -n "$used" ] || { echo "  $TAG: no seeds succeeded, not merging"; return; }
  echo "  $TAG: merging seeds $used"
  root -b -q -l "merge_acceptance.C(\"$TAG\",\"$used\",$EX,\"$OUT/\")" > "$OUT/${TAG}_merge.log" 2>&1
  mv -f "acceptance_merged_${TAG}.root" "$OUT/" 2>/dev/null || mv -f "$HERE/acceptance_merged_${TAG}.root" "$OUT/" 2>/dev/null || true
  rm -f "acceptance_merged_${TAG}.png" "$HERE/acceptance_merged_${TAG}.png"
  [ -s "$OUT/acceptance_merged_${TAG}.root" ] && echo "  $TAG: wrote $OUT/acceptance_merged_${TAG}.root" \
                                              || echo "  $TAG: MERGE_FAILED"
}

echo "vertex-z slab ${ZMIN}-${ZMAX} mm, chi2Cut $CHI2 -> $OUT"
run_level gs  0.0   "1001 1002 1003 1004 1005"
run_level ex1 6.094 "1001 1002 1003 1004 1005"
run_level ex8 8.317 "1001 1002 1003"
echo "done"
