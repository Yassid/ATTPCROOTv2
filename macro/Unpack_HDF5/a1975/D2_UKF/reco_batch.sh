#!/usr/bin/env bash
# Batch reco of the first 10 a1975 D2 (d,p)17C physics runs -> /mnt/f/a1975/reco_d2/
# Resumable (reco_one.sh skips finished runs). Parallelism configurable.
#   Usage: reco_batch.sh [Nparallel]   (default 4; drop to 2 if F:/drvfs I/O stalls)
# NOTE: proton-target production found drvfs caps at ~2 concurrent big-h5 reads;
#       4 is attempted per user request -- watch throughput, throttle to 2 if needed.
HERE="$(cd "$(dirname "$0")" && pwd)"
NPAR="${1:-4}"
RUNS=(run_0016 run_0017 run_0018 run_0019 run_0020 run_0021 run_0022 run_0023 run_0026 run_0027)

echo "=== D2 reco batch: ${#RUNS[@]} runs, ${NPAR}-parallel, $(date) ==="
printf '%s\n' "${RUNS[@]}" | xargs -P "$NPAR" -I{} bash "$HERE/reco_one.sh" {}
echo "=== batch finished, $(date) ==="
ls -lh /mnt/f/a1975/reco_d2/*_reco.root 2>/dev/null
