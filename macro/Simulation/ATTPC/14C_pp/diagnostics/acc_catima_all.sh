#!/usr/bin/env bash
# All three levels x five seeds of the CATIMA acceptance, N-way parallel.
#   ./acc_catima_all.sh [npar] [nEvents]
NPAR=${1:-4}; NEV=${2:-8000}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
JOBS=""
for s in 1001 1002 1003 1004 1005; do
  JOBS="$JOBS gs:0.0:$s ex1:6.094:$s ex8:8.317:$s"
done
printf "%s\n" $JOBS | xargs -P "$NPAR" -I{} bash -c '
  IFS=: read -r tag ex seed <<< "{}"
  '"$HERE"'/acc_catima_C14.sh "$tag" "$ex" "$seed" '"$NEV"' 2>&1 | tail -1'
echo "[$(date +%H:%M:%S)] ALL ACCEPTANCE JOBS DONE"
