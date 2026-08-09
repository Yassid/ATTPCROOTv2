#!/usr/bin/env bash
# Turn finished GENFIT acceptance jobs for one level into the merged no-chi2 acceptance that
# the correction macros consume.
#
# acc_batch_genfit.sh stops at the per-seed acceptance. Two more steps are needed and were
# previously run by hand, which is why they are easy to forget:
#
#   acceptance_split_C14.C   per seed, splits the numerator by which cut removed each event and
#                            stores the counterfactual hNoChi2_<level> alongside the denominator
#   export_nochi2_acc.C      sums the seeds and writes acceptance_merged_<level>.root
#
# The no-chi2 variant is the one used for the excited states: the chi2/ndf cut removes short
# backward tracks preferentially and puts a hole in the acceptance that is a property of the cut
# rather than of the apparatus.
#
#   ./acc_finish_level.sh <tag> <resEx> [seedLo] [seedHi]
#   ./acc_finish_level.sh ex8 8.317 1001 1003
set -eo pipefail
TAG=$1; EX=$2; S0=${3:-1001}; S1=${4:-1005}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); SIM=$(cd "$HERE/.." && pwd)
REPO=$(cd "$SIM/../../../.." && pwd)
ACC=${ACC_OUT:-/mnt/f/a1954_C14_acc_gf}
OUTDIR=${NOCHI2_OUT:-/mnt/f/a1954_C14_acc_gf_nochi2}
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
cd "$SIM"
mkdir -p diagnostics/split "$OUTDIR"

NDONE=0
for ((s=S0; s<=S1; ++s)); do
  J="${TAG}_s${s}"
  if [ ! -f "$ACC/$J.marker" ]; then
    echo "  $J: no COMPLETED marker -- skipping"
    continue
  fi
  root -b -q -l "acceptance_split_C14.C(\"$ACC/${J}_sim.root\",\"$ACC/${J}_genfit.root\",\"$TAG\",$EX,161.0,5.0,36,180.0,10.0,0.5,2.0,kTRUE,\"gf_s${s}\")" \
       > "$ACC/${J}_split.log" 2>&1
  if [ -s "diagnostics/split/acc_split_${TAG}_gf_s${s}.root" ]; then
    echo "  $J: split done"; NDONE=$((NDONE+1))
  else
    echo "  $J: SPLIT_FAILED -- see $ACC/${J}_split.log"
  fi
done

# Averaging over "however many happened to work" silently understates the seed spread, so refuse
# rather than produce a merged file whose seed count is not what was asked for.
[ "$NDONE" -gt 0 ] || { echo "no seeds split -- nothing to merge"; exit 1; }
echo "  merging $NDONE seed(s)"
root -b -q -l "export_nochi2_acc.C(\"diagnostics/split\",\"$TAG\",\"gf\",\"$OUTDIR\",$S0,$S1)"
echo "  wrote $OUTDIR/acceptance_merged_${TAG}.root"
