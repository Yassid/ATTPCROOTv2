#!/usr/bin/env bash
# Full D2-target 16C(d,p)17C production over the 47 canonical (d,p) runs.
# Pipelines reco (2-parallel, drvfs I/O-bound) with a genfit fit sweep (FITPAR-parallel,
# CPU-bound, leak-fixed + seed-fix on) that picks up runs as their reco completes.
# Both stages are resumable via .done markers, so this can be re-launched safely.
#   Usage: prod_full.sh [FITPAR]   (default 3)
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
FITPAR="${1:-3}"
RECOPAR=2
RECODIR=/mnt/f/a1975/reco_d2

# 47 canonical (d,p) runs (union of C16_pd_ana.C + _ang_dist.C; from C16_pd_AngBins.C)
NUMS="0016 0017 0018 0019 0020 0021 0022 0023 0026 0027 0031 0032 0034 0036 0037 0038 \
0039 0040 0041 0042 0043 0044 0046 0048 0057 0058 0076 0077 0078 0079 0080 0082 0083 \
0084 0085 0086 0087 0088 0089 0091 0092 0095 0096 0097 0098 0102 0103"
RUNS=""; for n in $NUMS; do RUNS="$RUNS run_$n"; done
NTOT=$(echo $RUNS | wc -w)

echo "=== D2 (d,p) FULL PRODUCTION: $NTOT runs, reco ${RECOPAR}-par + fit ${FITPAR}-par, $(date) ==="

# 1) reco sweep in the background (reco_one.sh skips runs already done)
printf '%s\n' $RUNS | xargs -P "$RECOPAR" -I{} bash "$HERE/reco_one.sh" {} &
RECOPID=$!
echo "[reco] sweep launched PID $RECOPID"

# 2) fit sweep loop: each pass fits every reco'd-but-unfit run (fit_one.sh skips done /
#    missing-reco), then waits; exits when all $NTOT fits are done, or reco ended + a
#    final mop-up pass is complete.
while true; do
  printf '%s\n' $RUNS | xargs -P "$FITPAR" -I{} bash "$HERE/fit_one.sh" {} -2.85
  nfit=$(ls $RECODIR/*_genfitter_p.root.done 2>/dev/null | wc -l)
  nreco=$(ls $RECODIR/*_reco.root.done 2>/dev/null | wc -l)
  echo "[sweep] $(date '+%H:%M:%S')  reco-done ${nreco}/${NTOT}  fit-done ${nfit}/${NTOT}"
  [ "$nfit" -ge "$NTOT" ] && { echo "=== ALL $NTOT FITS DONE $(date) ==="; break; }
  if ! kill -0 "$RECOPID" 2>/dev/null; then
    printf '%s\n' $RUNS | xargs -P "$FITPAR" -I{} bash "$HERE/fit_one.sh" {} -2.85
    nfit=$(ls $RECODIR/*_genfitter_p.root.done 2>/dev/null | wc -l)
    echo "=== reco ended; final fit pass done: ${nfit}/${NTOT} $(date) ==="
    break
  fi
  sleep 150
done
echo "=== PRODUCTION DRIVER EXIT $(date) ==="
