#!/usr/bin/env bash
# Per-run worker: full unpack+PSA+SC+PRA of one a1975 D2 run -> reco_d2/<run>_reco.root
# Resumable: skips if the output already exists and is >50 MB.
# Usage: reco_one.sh run_0016
set -e
RUN="$1"
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
OUTDIR=/mnt/f/a1975/reco_d2/
OUT="${OUTDIR}${RUN}_reco.root"

# skip if already COMPLETE. A full run's _reco.root is ~1GB; a partial (interrupted)
# file can still exceed tens of MB, so size alone is NOT a completion signal. Trust a
# .done marker (written only after root exits 0), with a >700MB size fallback for runs
# completed before this marker scheme existed.
if [ -f "${OUT}.done" ]; then echo "[skip] $RUN (done marker)"; exit 0; fi
if [ -f "$OUT" ]; then
  sz=$(stat -c %s "$OUT" 2>/dev/null || echo 0)
  if [ "$sz" -gt 734003200 ]; then echo "[skip] $RUN (>700MB, complete)"; exit 0; fi
fi

source "$REPO/build/config.sh" >/dev/null 2>&1
cd "$REPO/macro/Unpack_HDF5/a1975/D2_UKF"
echo "[start] $RUN  $(date '+%H:%M:%S')"
# nEvents=-1 (all), persistRaw=false (pattern-only, ~1GB), explicit outDir
root -l -b -q "unpackReco_a1975_deuterium.C(\"$RUN\", -1, false, \"$OUTDIR\")" \
   > "$REPO/macro/Unpack_HDF5/a1975/D2_UKF/log_${RUN}.txt" 2>&1
touch "${OUT}.done"   # reached only if root exited 0 (set -e aborts otherwise)
echo "[done ] $RUN  $(date '+%H:%M:%S')  -> $(du -h "$OUT" 2>/dev/null | cut -f1)"
