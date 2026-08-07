#!/usr/bin/env bash
# 5 seeds x 2 levels = 10 acceptance jobs, 3 at a time (8 cores; reco is single-threaded and
# drvfs writes to F: are not fast, so leave headroom rather than saturating).
set -eo pipefail
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
NEV=${1:-8000}
JOBS=""
for S in 1001 1002 1003 1004 1005; do
  JOBS+="gs 0.0 $S\n"
  JOBS+="ex1 6.094 $S\n"
done
printf "$JOBS" | xargs -P 3 -L 1 bash -c "$HERE/acc_batch.sh \$0 \$1 \$2 $NEV"
echo ACC_PARALLEL_COMPLETED
