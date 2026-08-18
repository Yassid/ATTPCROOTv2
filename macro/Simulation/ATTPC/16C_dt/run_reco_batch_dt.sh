#!/usr/bin/env bash
# Digitise + reconstruct all five 16C(d,t)15C states, at the gain matched to the data.
#
# GAIN 35000, AND IT IS MEASURED, NOT INHERITED. Gain decides how many electrons reach a pad,
# hence whether the pad crosses the PSA threshold, hence how many hits a track keeps -- which is
# exactly what the pattern finder and the PID see, and therefore what an acceptance measures.
# Scanned on gs_s3001 (200 events) against run_0016, non-beam tracks (>500 hits removed):
#
#   gain     median hits/track   sim/data
#   10000            39            0.51     <- the production par's value: HALF the data
#   25000            62            0.82
#   35000            75            0.99     <- adopted
#   50000            90            1.18
#   150000          118            1.55     <- the (p,d) simulation par's value
#   400000          135            1.78
#
# DATA run_0016 median = 76.
#
# DO NOT COPY THE (p,d) VALUE. Its simulation par declares 150000, and measured on the same
# footing that simulation gives median 117 hits/track against 46 in run_0106 -- a factor 2.54.
# That is the long-standing "sim/data diverge structurally" item, and it is gain, not clustering.
# Simulated tracks with 2.5x the hits are easier to find than real ones, so an acceptance measured
# that way is too high, and a yield divided by it comes out too low. Taking 150000 as a sensible
# default here would have reproduced the same inflation at 1.55x.
#
# Event efficiency saturates above ~25000 (99 of 200 events reconstruct, against 92 at 10000 and
# 100 at 150000), so gain is not choosing how many events survive -- it is choosing how much of
# each track survives. That is why the ratio, not the efficiency, is the thing to match.
#
# The scan is crude on purpose: all non-beam tracks on one run, not PID-selected tritons. It is
# strong enough to reject 10000 and 150000, NOT to defend 35000 to two digits. Re-measure against
# a PID-gated sample once the gate exists.
#
#   ./run_reco_batch_dt.sh [nparallel] [gain] [simDir]
set -uo pipefail
NPAR="${1:-3}"
GAIN="${2:-35000}"
SIM="${3:-/mnt/f/a1975_C16_dt_sim}"
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../../../.." && pwd)
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
export PAR="ATTPC.a1975_deuterium_dv1104_g${GAIN}.par"
export SIM HERE
cd "$HERE"

[ -s "$REPO/parameters/$PAR" ] || { echo "ERROR: no $PAR -- run ./make_gain_pars.sh $GAIN"; exit 1; }

one() {
  tag="$1"
  in="${SIM}/${tag}_s3001_sim.root"
  out="${SIM}/${tag}_s3001_reco.root"
  # generation writes its marker only on success; without this check a reco could run on a file
  # that a killed generator left half written, and nothing downstream would look wrong
  [ -f "${SIM}/${tag}_s3001.marker" ] || { echo "[nogen] $tag -- generation not marked complete"; return 0; }
  [ -s "$in" ] || { echo "[noinput] $tag"; return 0; }
  [ -f "${out}.done" ] && { echo "[have] $tag"; return 0; }
  root -b -l -q "run_reco_C16dt.C(\"$in\",\"$out\",\"$PAR\")" > "${SIM}/${tag}_reco.log" 2>&1
  if [ ! -s "$out" ]; then
    echo "[FAIL] $tag  (see ${SIM}/${tag}_reco.log)"; return 0
  fi
  # a run whose par silently did not load would reconstruct at the wrong gain and look fine
  grep -q "par = $PAR" "${SIM}/${tag}_reco.log" \
    || echo "[WARN] $tag: reco did not report $PAR"
  touch "${out}.done"; echo "[ok] $tag  $(date '+%H:%M:%S')  $(du -h "$out" | cut -f1)"
}
export -f one

echo "=== (d,t) sim reco: 5 states, gain $GAIN, $NPAR parallel ==="
echo "=== par $PAR (drift velocity, TBEntrance and ZPadPlane are the data's) ==="
printf '%s\n' gs ex1 ex2 ex3 ex4 | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "=== done: $(ls "$SIM"/*_reco.root 2>/dev/null | wc -l) reco files in $SIM ==="
