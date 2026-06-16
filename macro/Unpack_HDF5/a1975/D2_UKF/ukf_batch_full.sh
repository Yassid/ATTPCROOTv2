#!/usr/bin/env bash
# Full UKF batch over the 47 (d,p) runs with the backward fix (clusterDirSeed + minClus4,
# via fitUKF_a1975_deuterium.C defaults). CPU-bound, resumable via _genfitter_p_UKF.root.done.
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
NPAR="${1:-4}"
NUMS="0016 0017 0018 0019 0020 0021 0022 0023 0026 0027 0031 0032 0034 0036 0037 0038 0039 0040 0041 0042 0043 0044 0046 0048 0057 0058 0076 0077 0078 0079 0080 0082 0083 0084 0085 0086 0087 0088 0089 0091 0092 0095 0096 0097 0098 0102 0103"
RUNS=""; for n in $NUMS; do RUNS="$RUNS run_$n"; done
echo "=== UKF FULL batch (fix on): 47 runs, ${NPAR}-par, $(date) ==="
printf '%s\n' $RUNS | xargs -P "$NPAR" -I{} bash "$HERE/fit_one_ukf.sh" {}
echo "=== UKF batch finished $(date): $(ls /mnt/f/a1975/reco_d2/*_genfitter_p_UKF.root.done 2>/dev/null|wc -l)/47 done ==="
