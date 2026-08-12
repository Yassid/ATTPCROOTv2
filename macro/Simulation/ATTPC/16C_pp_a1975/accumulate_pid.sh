#!/usr/bin/env bash
# Accumulate PID statistics for redrawing the proton gate.
#
# Only reco + an UNGATED fit are needed: AtSpyralPID works off the pattern track, and the fit is
# run purely because AtGenfitter is what writes the AtPIDEvent branch. The quality of the fit is
# irrelevant here, and gating it would defeat the purpose -- the whole point is to see where the
# band lies before deciding what to keep.
#
# One seed per job. Seeds must differ or the samples are byte-identical copies.
#
#   ./accumulate_pid.sh <seed> [nEvents] [outDir]
#
# RESUMABLE PER STAGE. The first attempt at this was killed by an out-of-memory shutdown of the
# whole WSL VM with the reco stage a quarter done, and an all-or-nothing script would have thrown
# away three complete generations to get back to where it was. Each stage is skipped only on
# evidence that it actually finished:
#   generation  -- the sim file exists AND holds the full nEvents. A truncated or 272-byte stub
#                  does not count, which is exactly what the killed seeds left behind.
#   reco / fit  -- output exists, is non-empty, and its log carries the completion line. A written
#                  file is not a finished job.
# The reco output is KEPT rather than deleted at the end: with the pad traces no longer persisted
# it is a fraction of its old size, and keeping it means a later change to the fit or the gate
# costs a refit instead of another full reconstruction.
set -eo pipefail
SEED=${1:?need a seed}; NEV=${2:-12000}; OUT=${3:-/mnt/f/a1975_C16_pp_pid}
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
mkdir -p "$OUT"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
J=s$SEED
[ -f "$OUT/$J.marker" ] && { echo "$J already done"; exit 0; }

# --- generation ------------------------------------------------------------------------------
# entries==NEV is the completion test; the generator alternates beam and reaction events, so this
# is 2x the reaction count and half these entries carry no primary proton.
# A missing or corrupt sim file makes ROOT exit non-zero, and that is the ordinary case here --
# it is precisely how a killed generation announces itself -- so this probe must not trip set -e.
nsim=$(root -b -q -l -e "TFile*f=TFile::Open(\"$OUT/${J}_sim.root\");TTree*t=(f&&!f->IsZombie())?(TTree*)f->Get(\"cbmsim\"):nullptr;printf(\"N %lld\\n\",t?t->GetEntries():-1);" 2>/dev/null | awk '/^N /{print $2}' || true)
nsim=${nsim:--1}
if [ "$nsim" -eq "$NEV" ]; then
  echo "[$(date +%H:%M:%S)] $J gen already complete ($nsim events)"
else
  echo "[$(date +%H:%M:%S)] $J generating ($NEV events, had ${nsim:--1})"
  cd "$REPO/macro/Simulation/ATTPC/16C_pp_a1975"
  root -b -q -l "C16_pp_a1975_sim.C($NEV,2.,178.,\"TGeant4\",-28.5,\"$OUT/${J}_sim.root\",0.0,$SEED)" > "$OUT/${J}_gen.log" 2>&1
  grep -q "RNG seed requested: $SEED" "$OUT/${J}_gen.log" || { echo "$J SEED_NOT_APPLIED"; exit 1; }
fi

# --- reconstruction --------------------------------------------------------------------------
if [ -s "$OUT/${J}_reco.root" ] && grep -q "sim reco done" "$OUT/${J}_reco.log" 2>/dev/null; then
  echo "[$(date +%H:%M:%S)] $J reco already complete"
else
  echo "[$(date +%H:%M:%S)] $J reco"
  cd "$REPO/macro/Simulation/ATTPC/16C_pp_a1975"
  root -b -q -l "run_reco_C16_TC.C(\"$OUT/${J}_sim.root\",\"$OUT/${J}_reco.root\",\"ATTPC.a1975_C16_sim.par\",20,0,20.0)" > "$OUT/${J}_reco.log" 2>&1
  grep -q "sim reco done" "$OUT/${J}_reco.log" || { echo "$J RECO_FAILED"; exit 1; }
fi

# --- PID pass ----------------------------------------------------------------------------------
# NO FIT. AtPIDTask binds AtPatternEvent, exactly as the fitter does, so the PID landscape this
# gate is drawn from does not depend on any track being fitted -- the fit was only ever the vehicle
# that happened to carry AtPIDTask along. pidPass_a1975.C runs the task on its own, which is about
# a minute per sample against an hour for a full genfit pass over the same events.
#
# The field is left at pidPass's default 2.85, which is also AtPIDTask's built-in default and thus
# what the data production used -- fitGenfitter never calls SetBField. Rigidities from these
# samples are therefore on the same scale as the data plane the gate gets overlaid on.
if [ -s "$OUT/${J}_pid.root" ] && [ -f "$OUT/${J}_pid.done" ]; then
  echo "[$(date +%H:%M:%S)] $J pid already complete"
else
  echo "[$(date +%H:%M:%S)] $J pid"
  cd "$REPO/macro/Unpack_HDF5/a1975/UKF"
  root -b -q -l "pipeline/pidPass_a1975.C(\"$J\",\"$OUT/\",\"$OUT/\",2.85)" > "$OUT/${J}_pid.log" 2>&1
  [ -s "$OUT/${J}_pid.root" ] || { echo "$J PID_FAILED"; exit 1; }
  touch "$OUT/${J}_pid.done"
fi

echo COMPLETED > "$OUT/$J.marker"
echo "[$(date +%H:%M:%S)] $J done"
