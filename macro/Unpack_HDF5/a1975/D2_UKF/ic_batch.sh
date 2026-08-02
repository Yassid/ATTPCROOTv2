#!/usr/bin/env bash
# Ion-chamber (beam-particle) extraction for the 47 a1975 D2-target runs -> /mnt/f/a1975/ic_d2/
# <run>_IC.root (~1 MB each, ~2.5 min/run). Resumable: a run with an existing output is skipped.
#   ./ic_batch.sh [nparallel]
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
NPAR="${1:-4}"
OUTDIR="${2:-/mnt/f/a1975/ic_d2/}"
mkdir -p "$OUTDIR"
NUMS="0016 0017 0018 0019 0020 0021 0022 0023 0026 0027 0031 0032 0034 0036 0037 0038 0039 0040 0041 0042 0043 0044 0046 0048 0057 0058 0076 0077 0078 0079 0080 0082 0083 0084 0085 0086 0087 0088 0089 0091 0092 0095 0096 0097 0098 0102 0103"
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1

one() {
  n="$1"; out="$2/run_${n}_IC.root"
  # The macro unlinks its dummy only after out.Close(), so "dummy gone" is the success marker.
  # Testing -s alone re-skipped the 266-byte stubs left by runs that died mid-loop.
  [ -s "$out" ] && [ ! -f "$2/run_${n}_unpackIC_dummy.root" ] && { echo "[skip] run_$n"; return; }
  rm -f "$out" "$2/run_${n}_unpackIC_dummy.root"
  root -b -l -q "unpackIC_d2.C(\"run_${n}\",-1,\"$2\")" > "iclog_${n}.txt" 2>&1
  echo "[done] run_$n $(date '+%H:%M:%S') $(du -h "$out" 2>/dev/null | cut -f1)"
}
export -f one
echo "=== IC batch: $(echo $NUMS | wc -w) runs, ${NPAR}-par, start $(date) ==="
printf '%s\n' $NUMS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {} "$OUTDIR"
echo "=== IC batch done $(date): $(ls "$OUTDIR"/*_IC.root 2>/dev/null | wc -l) files ==="
