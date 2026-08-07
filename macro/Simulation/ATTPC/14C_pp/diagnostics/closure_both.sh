#!/usr/bin/env bash
# Fit the SAME truth-gated sim tracks with BOTH field-sign conventions and print the closure.
#
# The sim's truth is exactly Ebeam = 161.00 MeV, so the g.s. centroid mu IS the bias: the correct
# sign is the one that lands mu at 0. This is a measurement, not an argument -- a wrong bFieldSign
# does not crash, it silently biases KE, so it has to be decided on the closure and not on
# reasoning about conventions.
#
#   ./closure_both.sh [NAME] [INDIR]        default negBg, diagnostics/negB/
set -eo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); SIM=$(cd "$HERE/.." && pwd)
REPO=$(cd "$SIM/../../../.." && pwd); UKF="$REPO/macro/Unpack_HDF5/a1954/UKF"
NAME=${1:-negBg}; IN=${2:-"$HERE/negB/"}
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
[ -s "${IN%/}/${NAME}_reco.root" ] || { echo "MISSING ${IN%/}/${NAME}_reco.root"; exit 1; }
cd "$SIM"

# Symlink the one reco under two names so the two fits cannot overwrite each other's output.
ln -sf "${IN%/}/${NAME}_reco.root" "${IN%/}/${NAME}M_reco.root"
ln -sf "${IN%/}/${NAME}_reco.root" "${IN%/}/${NAME}P_reco.root"

for CFG in "M:-1" "P:+1"; do
  SUF=${CFG%%:*}; S=${CFG##*:}
  echo "[$(date +%H:%M:%S)] === UKF fit, bFieldSign = $S  (${NAME}${SUF})"
  root -b -q -l "$UKF/pipeline/fitUKF_C14.C(\"${NAME}${SUF}\",-1,\"proton\",$S,2.85,3.553e-5,\"\",\"$IN\",0.5,0.1,1,10,\"$IN\")" 2>&1 | tail -2
  root -b -q -l "$UKF/pp/ex_C14.C(\"${NAME}${SUF}\",\"$IN\",161.0,1e9,\"_${NAME}${SUF}\",1.007825,14.003242,\"\",\"ukf\")" 2>&1 | grep "good track"
done
echo "CLOSURE_BOTH_COMPLETED"
