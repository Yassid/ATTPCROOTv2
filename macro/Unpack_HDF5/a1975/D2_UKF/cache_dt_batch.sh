#!/usr/bin/env bash
# Build the full (d,t) kinematics cache: one ex_dt_a1975.C pass per run (PID plane +
# vertex + chi2 + run number all kept), then hadd into a single ntuple.
#   ./cache_dt_batch.sh [nparallel] [outdir]
# Resumable: a run whose per-run cache already exists is skipped.
HERE="$(cd "$(dirname "$0")" && pwd)"; cd "$HERE"
NPAR="${1:-4}"
OUTDIR="${2:-/tmp/claude-1000/-home-yassid/16bbcfdc-9a7e-4868-b7ba-0efb390913e3/scratchpad/dtcache}"
mkdir -p "$OUTDIR"
NUMS="0016 0017 0018 0019 0020 0021 0022 0023 0026 0027 0031 0032 0034 0036 0037 0038 0039 0040 0041 0042 0043 0044 0046 0048 0057 0058 0076 0077 0078 0079 0080 0082 0083 0084 0085 0086 0087 0088 0089 0091 0092 0095 0096 0097 0098 0102 0103"
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1

one() {
  n="$1"; out="$2/run_${n}.root"
  [ -s "$out" ] && { echo "[skip] run_$n"; return; }
  root -b -l -q "ex_dt_a1975.C(\"run_${n}\",\"/mnt/f/a1975/reco_d2/\",\"\",\"pid/triton_d2.json\",180.0,10.0,0.0,90.0,10.0,90.0,\"${out}\",\"/dev/null.png\")" \
     > "$2/log_${n}.txt" 2>&1
  echo "[done] run_$n $(date '+%H:%M:%S')"
}
export -f one
printf '%s\n' $NUMS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {} "$OUTDIR"
hadd -f "$OUTDIR/../dt_kin_full.root" "$OUTDIR"/run_*.root > "$OUTDIR/hadd.log" 2>&1
echo "=== cache built: $OUTDIR/../dt_kin_full.root ==="
