#!/usr/bin/env bash
# C15d POST-RECO chain for the INTEGRAL-dE/dx reconstruction in
# /home/yassid/Data/C15d_reco_mfi (built by rereco_mfi_C15d.sh).
# Takes the chain as far as it can go WITHOUT a gate and then stops, per the standing rule that
# fitting before gating wasted ~14x the CPU last time.
#
#   ./postreco_hdb_C15d.sh [nparallel]
#
# Stages:
#   1  pidntuple_C15d.C per run      -> <run>_pid.root beside each reco  (a few MB each)
#      (its 3rd argument is the OUTPUT dir, not the IC dir -- the usage line in that macro
#       used to say ~/C15d_ic/ and the cache went there instead)
#   2  measure_gain_C15d.C           -> gainmatch_mfi_C15d.csv
#   3  pid/make_points_C15d.C        -> pid/points_mfi_C15d.root         (the gate GUI's plane)
#   -- STOP. Draw the gate, then run the fits with C15D_RECO pointed here.
#
# ★ STAGE 2 IS NOT OPTIONAL AND THE OLD TABLE MUST NOT BE REUSED. `gainmatch_C15d.csv` was
#   measured on the AtPSAMax reco. MultiFit resolves overlapping pulses, so the hits each track
#   carries change (126 -> 138 hits/track) and so does their charge; dE/dx is computed from those
#   hits, so the per-run gains are a different set of numbers. Drawing a gate on a plane carrying
#   the old table smears it by the per-run drift (~29 %).
#
# The IC summaries are NOT rebuilt: ic_batch.sh derives them from the raw /frib group, which the
# reconstruction never touches, so the 103 existing <run>_ic.root files stay valid.
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
RECO="${C15D_RECO:-/home/yassid/Data/C15d_reco_mfi}"
IC="${C15D_IC:-/home/yassid/C15d_ic}"
NPAR="${1:-10}"
GAINCSV=gainmatch_mfi_C15d.csv
POINTS=pid/points_mfi_C15d.root
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
cd "$HERE"

mapfile -t RUNS < <(ls -1 "$RECO"/run_*_reco.root 2>/dev/null | sed -E 's/.*run_0*([0-9]+)_reco\.root/\1/')
echo "=== post-reco chain on $RECO : ${#RUNS[@]} reconstructed runs ==="
[ ${#RUNS[@]} -gt 0 ] || { echo "nothing to do"; exit 1; }

# ---- stage 1: per-run PID caches
pidone(){ local i="$1"; local r; r=$(printf run_%04d "$i")
  [ -s "$RECO/${r}_pid.root" ] && return 0
  root -b -q "$HERE/pidntuple_C15d.C(\"$r\",\"$RECO/\",\"$RECO/\")" >"$RECO/logs/${r}_pid.log" 2>&1 || true
  [ -s "$RECO/${r}_pid.root" ] && echo "  $r pid cache ok" || echo "  $r PID CACHE FAILED"
}
export -f pidone; export RECO IC HERE
mkdir -p "$RECO/logs"
echo "--- stage 1: PID caches ($NPAR parallel)"
printf '%s\n' "${RUNS[@]}" | xargs -P "$NPAR" -I{} bash -c 'pidone "$@"' _ {}
NC=$(ls -1 "$RECO"/*_pid.root 2>/dev/null | wc -l)
echo "--- $NC PID caches present"

# ---- stage 2: re-measure the per-run gain on THIS reco
echo "--- stage 2: gain match -> $GAINCSV"
root -b -q "measure_gain_C15d.C(0.30,0.60,\"$RECO/\",\"$GAINCSV\",200,200,0,0,0,0,-1,\"plots_mfi/\")" 2>&1 \
  | grep -vE '^Warning|^Processing|libGL' | tail -25
[ -s "$GAINCSV" ] || { echo "ERROR: $GAINCSV not written"; exit 1; }

# ---- stage 3: the points file the gate GUI reads
echo "--- stage 3: points file -> $POINTS"
root -b -q "pid/make_points_C15d.C(\"$RECO/\",\"$POINTS\",\"$GAINCSV\",\"$IC/\")" 2>&1 \
  | grep -vE '^Warning|^Processing|libGL' | tail -30
[ -s "$POINTS" ] || { echo "ERROR: $POINTS not written"; exit 1; }

cat <<MSG

=== STOPPING HERE ON PURPOSE ===
The plane is ready:  $POINTS
Draw the deuteron gate:

    cd $HERE && ./pid/open_gate_draw.sh deuteron $POINTS

then fit ONLY the gated tracks:

    C15D_RECO=$RECO C15D_FIT=/home/yassid/Data/C15d_fit_mfi \\
      ./fit_batch.sh <nparallel> pid/deuteron_C15d.json

⚠ open_gate_draw.sh DEFAULTS to pid/points_C15d.root (the OLD TriplClust plane). Pass the path above, or the
  gate will be drawn on the AtPSAMax plane and applied to the MultiFit one.
MSG
