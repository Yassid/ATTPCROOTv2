#!/usr/bin/env bash
# GENFIT half of the both-signs closure (closure_both.sh does UKF). Independent fitter, same
# tracks, same two conventions -- the two fitters agreeing is what makes the answer trustworthy.
# GENFIT takes the SIGNED field directly (not a +-1 flag): data = -2.85, legacy sim = +2.85.
set -eo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); SIM=$(cd "$HERE/.." && pwd)
REPO=$(cd "$SIM/../../../.." && pwd); UKF="$REPO/macro/Unpack_HDF5/a1954/UKF"
NAME=${1:-negBg}; IN=${2:-"$HERE/negB/"}
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
cd "$SIM"
for CFG in "M:-2.85" "P:+2.85"; do
  SUF=${CFG%%:*}; B=${CFG##*:}
  [ -e "${IN%/}/${NAME}${SUF}_reco.root" ] || ln -sf "${IN%/}/${NAME}_reco.root" "${IN%/}/${NAME}${SUF}_reco.root"
  echo "[$(date +%H:%M:%S)] === GENFIT, bField = $B  (${NAME}${SUF})"
  root -b -q -l "$UKF/pipeline/fitGenfit_C14.C(\"${NAME}${SUF}\",-1,\"$IN\",\"\",\"$IN\",$B,2,5,\"\",4.0,10.0,170.0,kFALSE,kFALSE,\"proton\",\"ATTPC_H300torr\")" 2>&1 | tail -2
  root -b -q -l "$UKF/pp/ex_C14.C(\"${NAME}${SUF}\",\"$IN\",161.0,1e9,\"_${NAME}${SUF}_gf\",1.007825,14.003242,\"\",\"genfit\")" 2>&1 | grep "good track"
done
echo "CLOSURE_GENFIT_COMPLETED"
