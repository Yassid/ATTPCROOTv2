#!/usr/bin/env bash
# Full UKF batch, TRITON hypothesis, over the 47 D2 runs -> _genfitter_t_UKF.root.
# Resumable: completed runs carry a .done marker and are skipped; a .root without a
# .done is treated as an interrupted leftover and re-run from scratch.
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
NPAR="${1:-4}"
NUMS="0016 0017 0018 0019 0020 0021 0022 0023 0026 0027 0031 0032 0034 0036 0037 0038 0039 0040 0041 0042 0043 0044 0046 0048 0057 0058 0076 0077 0078 0079 0080 0082 0083 0084 0085 0086 0087 0088 0089 0091 0092 0095 0096 0097 0098 0102 0103"
NTOT=$(echo $NUMS | wc -w)
RUNS=""; for n in $NUMS; do RUNS="$RUNS run_$n"; done
echo "=== UKF-t batch: $NTOT runs, ${NPAR}-par, start $(date) ==="
echo "=== already done: $(ls /mnt/f/a1975/reco_d2/*_genfitter_t_UKF.root.done 2>/dev/null|wc -l)/$NTOT ==="
printf '%s\n' $RUNS | xargs -P "$NPAR" -I{} bash "$HERE/fit_one_ukf_t.sh" {}
echo "=== UKF-t batch finished $(date): $(ls /mnt/f/a1975/reco_d2/*_genfitter_t_UKF.root.done 2>/dev/null|wc -l)/$NTOT done ==="
