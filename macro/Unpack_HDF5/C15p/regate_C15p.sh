#!/usr/bin/env bash
# Re-run everything downstream of a redrawn gate, and report whether it actually improved.
#
#   ./regate_C15p.sh [gate.json] [species] [nparallel]
#
# Gate -> fits -> kinematics -> Ex fit. Prints the enrichment (is it the right band?) and the
# g.s. / cluster widths (did the contamination go?), so a gate change can be judged on the
# spectrum rather than on how the polygon looks.
#
# The g.s. is the ONLY isolated level in this spectrum, so it is the only place a width means
# resolution. The 6-7 MeV structure is the unresolved 14C cluster (6.09/6.59/6.73/7.01) and its
# width is an envelope, not a resolution -- do not read it as one.
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
GATE="${1:-pid/deuteron_15Cp.json}"
SP="${2:-d}"
NPAR="${3:-20}"
RECO=/home/yassid/Data/a2091_C15p_reco_hdb
FIT=/home/yassid/Data/a2091_C15p_fit_hdb
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
cd "$HERE"
[[ -s "$GATE" ]] || { echo "ERROR: no gate at $GATE" >&2; exit 1; }

echo "=== 1/3  clearing the previous $SP fits (a gate change invalidates every one) ==="
rm -f "$FIT"/*_genfit_${SP}.root "$FIT"/*_kin_${SP}.root "$FIT"/in_${SP}"/"*.root 2>/dev/null || true

echo "=== 2/3  gated fits ==="
C15P_RECO=$RECO C15P_FIT=$FIT C15P_LOGS=/home/yassid/C15p_hdb_logs \
  ./fit_batch.sh "$NPAR" "$SP" "$HERE/runs_pp_hdb.txt" -1 "$GATE" 2>&1 | tail -3

echo "=== 3/3  kinematics + Ex ==="
root -b -q "pp/kin_pp_C15p.C(\"$FIT/\",\"\",\"pp/plots_hdb/\",5.0,60.0,\"pid/points_hdb_C15p.root\",\"$SP\")" 2>&1 \
  | grep -E "fitted|in gate|backward" || true
root -b -q "pp/exfit_C15p.C(\"pp/plots_hdb/pp_kin_C15p_${SP}.root\",198.0)" 2>&1 \
  | grep -E "entries|G.S.|CLUSTER|spacing" || true
echo
echo "compare against the previous gate before keeping it; v1 was G.S. FWHM 0.911, cluster 1.543, 18264 tracks."
