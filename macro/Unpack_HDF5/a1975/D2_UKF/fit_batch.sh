#!/usr/bin/env bash
# Batch proton-hypothesis genfit over the first 10 D2 reco files -> reco_d2/<run>_genfitter_p.root
# CPU-bound (reads small _reco.root) so 4-parallel is fine. Resumable.
#   Usage: fit_batch.sh [Nparallel] [bField]   (defaults: 4, -2.85)
HERE="$(cd "$(dirname "$0")" && pwd)"
NPAR="${1:-4}"
BFIELD="${2:--2.85}"
RUNS=(run_0016 run_0017 run_0018 run_0019 run_0020 run_0021 run_0022 run_0023 run_0026 run_0027)

echo "=== D2 proton-hyp fit batch: ${#RUNS[@]} runs, ${NPAR}-parallel, B=${BFIELD}, $(date) ==="
printf '%s\n' "${RUNS[@]}" | xargs -P "$NPAR" -I{} bash "$HERE/fit_one.sh" {} "$BFIELD"
echo "=== fit batch finished, $(date) ==="
ls -lh /mnt/f/a1975/reco_d2/*_genfitter_p.root 2>/dev/null
