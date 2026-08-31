#!/usr/bin/env bash
# One 46Ar(3He,d)47K sample: generate -> digitise+reconstruct -> PID observables. NO FIT.
#
#   ./accumulate_3Hed.sh <state> <seed> [nEvents] [outDir] [bFieldT] [padSize_mm] [simDir]
#
# <state> is gs | 360 | 2020, i.e. the 47K level: 1/2+ ground state, 0.36 MeV 3/2+, 2.02 MeV 7/2-.
# EACH STATE IS ITS OWN JOB with its own seed, so the three never share a random sequence and any
# one of them can be regenerated without touching the others. The tag is <state>_s<seed>, e.g.
# gs_s3001, and every file of a sample carries it.
#
# THE CHAIN STOPS AT THE PID OBSERVABLES, DELIBERATELY. AtPIDTask computes brho and sqrt(dE/dx)
# per pattern track and writes them out; it selects nothing. The deuteron gate gets drawn by hand
# on that plane (pid_plane_3Hed.C) before any fitting happens, so no stage here decides a species.
#
# RESUMABLE PER STAGE, on evidence that a stage actually finished -- a written file is not a
# finished job. The a1975 version of this script was written after an out-of-memory shutdown of
# the whole WSL VM threw away three complete generations, and the same reasoning applies:
#   generation  -- the sim file exists AND holds the full nEvents. A truncated or 272-byte stub
#                  does not count, and that is exactly what a killed generation leaves behind.
#   reco / pid  -- output exists, is non-empty, and its log carries the completion line.
set -eo pipefail
STATE=${1:?need a state: gs | 360 | 2020}; SEED=${2:?need a seed}
NEV=${3:-12000}; OUT=${4:-/mnt/f/ar46_3hed}
BT=${5:-2.85}          # tesla; the generator wants kG and the DATA sign convention (negative)
PAD=${6:--1}           # mm; <=0 keeps the real AT-TPC pad plane
SIMDIR=${7:-$OUT}      # generation depends on the FIELD only, so a 2 mm run reuses the sims
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF

case "$STATE" in
   gs)   EX=0.0   ;;
   360)  EX=0.360 ;;
   2020) EX=2.020 ;;
   *) echo "unknown state '$STATE' (want gs, 360 or 2020)"; exit 2 ;;
esac

BKG=$(awk -v b="$BT" 'BEGIN{printf "%.1f", -10*b}')
# PAR IS SELECTED BY FIELD, AND AN UNKNOWN FIELD IS A HARD ERROR. Each field now has its own
# par because DriftVelocity and CoefL/CoefT come from a Magboltz 300 torr scan AT THAT FIELD
# (make_3Hed_pars.sh) rather than being the H2 placeholders shared by every arm. A ?: fallback
# to the 2.85 T file would silently run a 2.0 T sim on 2.85 T transport and present as valid.
case "$(awk -v b="$BT" 'BEGIN{printf "%.2f", b}')" in
   2.00) PAR=ATTPC.46Ar_3Hed_sim_B20.par  ;;
   2.85) PAR=ATTPC.46Ar_3Hed_sim_B285.par ;;
   3.80) PAR=ATTPC.46Ar_3Hed_sim_B38.par  ;;
   *) echo "no par for B = $BT T -- add it to make_3Hed_pars.sh and re-run that script"; exit 2 ;;
esac
[ -f "$REPO/parameters/$PAR" ] || { echo "MISSING PAR $PAR -- run make_3Hed_pars.sh"; exit 2; }
mkdir -p "$OUT" "$SIMDIR"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
J=${STATE}_s${SEED}
echo "[cfg] B = $BT T ($BKG kG), pads = $PAD mm, par = $PAR, sim from $SIMDIR, out to $OUT"
[ -f "$OUT/$J.marker" ] && { echo "$J already done"; exit 0; }
cd "$REPO/macro/Simulation/ATTPC/46Ar_3Hed"

# --- generation ------------------------------------------------------------------------------
# entries == NEV is the completion test. The generator alternates beam and reaction events, so
# NEV entries is NEV/2 reactions and half the entries carry no deuteron -- do not read that as a
# 50 % efficiency. A missing or corrupt sim file makes ROOT exit non-zero, which is the ordinary
# way a killed generation announces itself, so this probe must not trip set -e.
nsim=$(root -b -q -l -e "TFile*f=TFile::Open(\"$SIMDIR/${J}_sim.root\");TTree*t=(f&&!f->IsZombie())?(TTree*)f->Get(\"cbmsim\"):nullptr;printf(\"N %lld\\n\",t?t->GetEntries():-1);" 2>/dev/null | awk '/^N /{print $2}' || true)
nsim=${nsim:--1}
if [ "$nsim" -eq "$NEV" ]; then
  echo "[$(date +%H:%M:%S)] $J gen already complete ($nsim entries)"
else
  echo "[$(date +%H:%M:%S)] $J generating: Ex = $EX MeV, $NEV entries, seed $SEED (had ${nsim:--1})"
  root -b -q -l "Ar46_3Hed_sim.C($NEV,15.,80.,\"TGeant4\",$BKG,\"$SIMDIR/${J}_sim.root\",$EX,$SEED)" \
       > "$SIMDIR/${J}_gen.log" 2>&1
  grep -q "RNG seed requested: $SEED" "$SIMDIR/${J}_gen.log" || { echo "$J SEED_NOT_APPLIED"; exit 1; }
  # Compare the excitation NUMERICALLY. ROOT prints 0.36 for 0.360, so a string match on the
  # argument as written fails on a run that is perfectly correct -- which is exactly what it did
  # the first time this script was run.
  exlog=$(awk '/47K excitation =/{print $4; exit}' "$SIMDIR/${J}_gen.log")
  awk -v a="${exlog:-nan}" -v b="$EX" 'BEGIN{exit !(a==a+0 && (a-b<1e-6 && b-a<1e-6))}' \
     || { echo "$J WRONG_STATE (log says '${exlog:-none}', wanted $EX)"; exit 1; }
fi

# --- digitisation + reconstruction -------------------------------------------------------------
if [ -s "$OUT/${J}_reco.root" ] && grep -q "sim reco done" "$OUT/${J}_reco.log" 2>/dev/null; then
  echo "[$(date +%H:%M:%S)] $J reco already complete"
else
  echo "[$(date +%H:%M:%S)] $J reco"
  root -b -q -l "run_reco_Ar46_TC.C(\"$SIMDIR/${J}_sim.root\",\"$OUT/${J}_reco.root\",\"$PAR\",20,0,20.0,15.0,7.5,\"tc\",20,$PAD)" \
       > "$OUT/${J}_reco.log" 2>&1
  grep -q "sim reco done" "$OUT/${J}_reco.log" || { echo "$J RECO_FAILED"; exit 1; }
fi

# --- PID observables (no gate, no fit) ---------------------------------------------------------
if [ -s "$OUT/${J}_pid.root" ] && grep -q "pid pass done" "$OUT/${J}_pid.log" 2>/dev/null; then
  echo "[$(date +%H:%M:%S)] $J pid already complete"
else
  echo "[$(date +%H:%M:%S)] $J pid"
  root -b -q -l "pidPass_Ar46.C(\"$J\",\"$OUT/\",\"$OUT/\",$BT,\"$PAR\")" > "$OUT/${J}_pid.log" 2>&1
  grep -q "pid pass done" "$OUT/${J}_pid.log" || { echo "$J PID_FAILED"; exit 1; }
fi

echo COMPLETED > "$OUT/$J.marker"
echo "[$(date +%H:%M:%S)] $J done"
