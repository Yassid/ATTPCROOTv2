#!/usr/bin/env bash
# Per-run worker: proton-hypothesis genfit of one D2 reco -> reco_d2/<run>_genfitter_p.root
# Reads the (smaller, ~1GB) <run>_reco.root, so this is CPU-bound -> safe at 4-parallel.
# Resumable: skips if the output already exists and is >5 MB.
# Usage: fit_one.sh run_0016 [bField]
set -e
RUN="$1"
BFIELD="${2:--2.85}"
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
IODIR=/mnt/f/a1975/reco_d2/
OUT="${IODIR}${RUN}_genfitter_p.root"

if [ ! -f "${IODIR}${RUN}_reco.root" ]; then
  echo "[miss] $RUN has no reco yet (${IODIR}${RUN}_reco.root) — skipping"
  exit 0
fi
# trust a .done marker (written only after root exits 0); a partial fit output can
# also exceed a few MB, so size alone is not a completion signal.
if [ -f "${OUT}.done" ]; then echo "[skip] $RUN fit (done marker)"; exit 0; fi

source "$REPO/build/config.sh" >/dev/null 2>&1
export LD_LIBRARY_PATH="$GENFIT/lib:$LD_LIBRARY_PATH"
cd "$REPO/macro/Unpack_HDF5/a1975/D2_UKF"
echo "[start] fit $RUN B=$BFIELD  $(date '+%H:%M:%S')"
root -l -b -q "fitGenfitter_a1975_deuterium.C(\"$RUN\", -1, \"$IODIR\", \"\", \"$IODIR\", $BFIELD)" \
   > "fitlog_${RUN}.txt" 2>&1
touch "${OUT}.done"   # reached only if root exited 0 (set -e aborts otherwise)
echo "[done ] fit $RUN  $(date '+%H:%M:%S')  -> $(du -h "$OUT" 2>/dev/null | cut -f1)"
