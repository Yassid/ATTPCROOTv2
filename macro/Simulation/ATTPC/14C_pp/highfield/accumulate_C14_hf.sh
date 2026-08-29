#!/usr/bin/env bash
# ONE SAMPLE of the 14C(p,p') field x pad-pitch study: generate -> reco -> genfit -> acceptance
# + Ex resolution, for a single (level, field, pitch, seed).
#
#   ./accumulate_C14_hf.sh <state> <B_T> <pad_mm> <seed> [nEvents]
#   ./accumulate_C14_hf.sh gs 7.0 2.0 7021
#
# <state> is one of the tags in the case block below: gs | ex6094 | ex6728 | ex7012 | ex8317.
#
# THE SIMS ARE SHARED PER FIELD, NOT PER CONFIGURATION. Transport depends on the magnetic field
# but not at all on the pad plane, so the AT-TPC-pitch and 2 mm runs of the same field read the
# SAME <state>_s<seed>_sim.root. A difference between them is then the pad plane and nothing
# else -- the same argument that made the 46Ar(3He,d) matrix interpretable. It also halves the
# generation cost. Two jobs of the same field can therefore race on the sim file, so generation
# takes a lock and the other job waits for it rather than writing over it.
#
# EVERYTHING ELSE IS THE ADOPTED 2026-08-25 a1954 RECIPE, held fixed: RT gas (3.308e-5), genfit
# with matEffects + native CATIMA, matFallback off, no manual gap eloss, acceptance at
# chi2/ndf < 5 on GetKinematicsXtr. The only things that vary across the campaign are the field,
# the pitch, and the diffusion coefficients that follow from the field (see the par files).
#
# RESUMABLE PER STAGE, on evidence a stage finished rather than on a file existing -- a killed
# job leaves a file behind. Stage tests: the sim must hold the full nEvents; every later stage
# must have written the completion line its log is checked for.
set -eo pipefail

STATE=${1:?need a state: gs | ex6094 | ex6728 | ex7012 | ex8317}
BT=${2:?need the field in tesla, e.g. 2.85}
PAD=${3:?need the pad pitch in mm, or -1 for the real AT-TPC pad plane}
SEED=${4:?need a seed}
NEV=${5:-8000}

case "$STATE" in
   gs)     EX=0.0   ;;   # elastic
   ex6094) EX=6.094 ;;   # 1-
   ex6728) EX=6.728 ;;   # 3-
   ex7012) EX=7.012 ;;   # 2+, the state that carries the B(E2)
   ex8317) EX=8.317 ;;   # 2+, isolated above the multiplet
   *) echo "unknown state '$STATE'"; exit 2 ;;
esac

REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
SIM="$REPO/macro/Simulation/ATTPC/14C_pp"
HF="$SIM/highfield"
UKF="$REPO/macro/Unpack_HDF5/a1954/UKF"
ROOTDIR=${HF_ROOT:-/mnt/f/a1954_C14_hf}

# tags: field to a 3-digit centi-tesla code, pitch to a name. b285_attpc, b700_2mm, ...
BTAG=$(awk -v b="$BT" 'BEGIN{printf "b%03d", b*100}')
PTAG=$(awk -v p="$PAD" 'BEGIN{ if (p>0) printf "%gmm", p; else printf "attpc" }')
CFG="${BTAG}_${PTAG}"
SIMDIR="$ROOTDIR/sims_$BTAG"
OUT="$ROOTDIR/$CFG"
mkdir -p "$SIMDIR" "$OUT"

# The par carries BField and the field-dependent diffusion coefficients. DIFFMODE=fixed selects
# the control set that holds CoefL/CoefT at the legacy 9e-4 for every field, which is what
# isolates "the field bent the tracks better" from "the field stopped the electrons spreading".
# 2.85 T IS THE ANCHOR ITSELF, so its "fixed diffusion" par would be its ordinary one; there is
# no _fixdiff file for it and asking for one is not an error.
PARSUF=""
[ "${DIFFMODE:-mb}" = "fixed" ] && [ "$BTAG" != "b285" ] && PARSUF="_fixdiff"
PAR="ATTPC.a1954_C14_hf_${BTAG}${PARSUF}.par"
[ -f "$REPO/parameters/$PAR" ] || { echo "MISSING PAR $PAR"; exit 2; }

BKG=$(awk -v b="$BT" 'BEGIN{printf "%.1f", -10*b}')   # kG, DATA sign convention (negative)
BNEG=$(awk -v b="$BT" 'BEGIN{printf "%.2f", -b}')     # tesla, for the fitter
J="${STATE}_s${SEED}"
JC="${J}_${CFG}"        # per-configuration tag: the reco/fit/acceptance products carry it

set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"

echo "[cfg] $JC : Ex = $EX MeV, B = $BT T ($BKG kG), pads = $PAD mm, par = $PAR, $NEV events"
[ -f "$OUT/$JC.marker" ] && { echo "$JC already COMPLETED"; exit 0; }

# ---- generation (shared per field) -------------------------------------------------------------
# Run it in a scratch cwd of its own: C14_pp_sim.C hardcodes ./data/attpcpar.root and
# ./data/geofile_C14_pp_full.root, so two generations sharing a working directory overwrite each
# other's parameter file.
nsim=$(root -b -q -l -e "TFile*f=TFile::Open(\"$SIMDIR/${J}_sim.root\");TTree*t=(f&&!f->IsZombie())?(TTree*)f->Get(\"cbmsim\"):nullptr;printf(\"N %lld\\n\",t?t->GetEntries():-1);" 2>/dev/null | awk '/^N /{print $2}' || true)
nsim=${nsim:--1}
if [ "$nsim" -ne "$NEV" ]; then
  LOCK="$SIMDIR/${J}.lock"
  if mkdir "$LOCK" 2>/dev/null; then
    trap 'rmdir "$LOCK" 2>/dev/null || true' EXIT
    WORK="$SIMDIR/work_$J"; mkdir -p "$WORK/data"
    echo "[$(date +%H:%M:%S)] $J generating (had $nsim)"
    ( cd "$WORK" && root -b -q -l "$SIM/C14_pp_sim.C($NEV,2.0,178.0,\"TGeant4\",$BKG,\"$SIMDIR/${J}_sim.root\",$EX,$SEED,\"ATTPC_H300torr_RT.root\")" ) > "$SIMDIR/${J}_gen.log" 2>&1
    grep -q "ATTPC_H300torr_RT" "$SIMDIR/${J}_gen.log" || { echo "$J WRONG_GAS"; exit 1; }
    grep -q "RNG seed requested: $SEED" "$SIMDIR/${J}_gen.log" || { echo "$J SEED_NOT_APPLIED"; exit 1; }
    # the field really has to have reached the transport, or every arm of the matrix is 2.85 T
    grep -q -- "$BKG" "$SIMDIR/${J}_gen.log" || echo "$J WARN: could not confirm $BKG kG in the log"
    rm -rf "$WORK"
    rmdir "$LOCK" 2>/dev/null || true; trap - EXIT
  else
    echo "[$(date +%H:%M:%S)] $J generation held by another job, waiting"
    for _ in $(seq 1 720); do sleep 10; [ -d "$LOCK" ] || break; done
    [ -d "$LOCK" ] && { echo "$J GEN_LOCK_TIMEOUT"; exit 1; }
  fi
  nsim=$(root -b -q -l -e "TFile*f=TFile::Open(\"$SIMDIR/${J}_sim.root\");TTree*t=(f&&!f->IsZombie())?(TTree*)f->Get(\"cbmsim\"):nullptr;printf(\"N %lld\\n\",t?t->GetEntries():-1);" 2>/dev/null | awk '/^N /{print $2}' || true)
  [ "${nsim:--1}" -eq "$NEV" ] || { echo "$J GEN_INCOMPLETE ($nsim of $NEV)"; exit 1; }
else
  echo "[$(date +%H:%M:%S)] $J generation already complete ($nsim entries)"
fi

cd "$SIM"

# ---- digitisation + reconstruction (per configuration) -----------------------------------------
# persistRaw = kFALSE: the pad traces are 4.7 GB per 8000 events and nothing past PSA reads them.
if [ -s "$OUT/${JC}_reco.root" ] && grep -q "sim reco done" "$OUT/${JC}_reco.log" 2>/dev/null; then
  echo "[$(date +%H:%M:%S)] $JC reco already complete"
else
  echo "[$(date +%H:%M:%S)] $JC reco"
  root -b -q -l "run_reco_C14.C(\"$SIMDIR/${J}_sim.root\",\"$OUT/${JC}_reco.root\",\"$PAR\",20,20,8,0,30.0,\"mover\",$PAD,500.0,kFALSE)" > "$OUT/${JC}_reco.log" 2>&1
  grep -q "sim reco done" "$OUT/${JC}_reco.log" || { echo "$JC RECO_FAILED"; exit 1; }
  if [ "$(awk -v p="$PAD" 'BEGIN{print (p>0)?1:0}')" = "1" ]; then
    grep -q "PADS  : square" "$OUT/${JC}_reco.log" || { echo "$JC PAD_PLANE_NOT_APPLIED"; exit 1; }
  fi
fi

# ---- genfit : matEffects + native CATIMA, the adopted configuration ------------------------------
if [ -s "$OUT/${JC}_genfit.root" ] && grep -q "dE/dx from CATIMA" "$OUT/${JC}_fit.log" 2>/dev/null; then
  echo "[$(date +%H:%M:%S)] $JC fit already complete"
else
  echo "[$(date +%H:%M:%S)] $JC genfit at B = $BNEG T"
  root -b -q -l "$UKF/pipeline/fitGenfit_C14.C(\"$JC\",-1,\"$OUT/\",\"\",\"$OUT/\",$BNEG,2,5,\"\",4.0,10.0,170.0,kTRUE,kFALSE,\"proton\",\"ATTPC_H300torr_RT\",kFALSE,kTRUE,0.0,1,kTRUE,kFALSE,kFALSE,kFALSE)" > "$OUT/${JC}_fit.log" 2>&1
  [ -s "$OUT/${JC}_genfit.root" ] || { echo "$JC FIT_FAILED"; exit 1; }
  grep -q "dE/dx from CATIMA" "$OUT/${JC}_fit.log" || { echo "$JC CATIMA_NOT_ENABLED"; exit 1; }
fi

# ---- acceptance vs theta_cm ---------------------------------------------------------------------
if ! grep -q "overall acceptance" "$OUT/${JC}_acc.log" 2>/dev/null; then
  echo "[$(date +%H:%M:%S)] $JC acceptance"
  root -b -q -l "$SIM/acceptance_C14.C(\"$SIMDIR/${J}_sim.root\",\"$OUT/${JC}_genfit.root\",\"$JC\",$EX,159.75,5.0,36,180.0,10.0,0.5,2.0,kTRUE)" > "$OUT/${JC}_acc.log" 2>&1
  grep -q "overall acceptance" "$OUT/${JC}_acc.log" || { echo "$JC ACC_FAILED"; exit 1; }
  mv -f "$SIM/diagnostics/acceptance_${JC}.root" "$OUT/" 2>/dev/null || true
  rm -f "$SIM/diagnostics/acceptance_${JC}.png"
fi

# ---- excitation-energy resolution ---------------------------------------------------------------
if ! grep -q "ex res done" "$OUT/${JC}_exres.log" 2>/dev/null; then
  echo "[$(date +%H:%M:%S)] $JC Ex resolution"
  root -b -q -l "$HF/ex_res_C14_hf.C(\"$SIMDIR/${J}_sim.root\",\"$OUT/${JC}_genfit.root\",\"$JC\",$EX,159.75,5.0,kTRUE,\"$OUT/\")" > "$OUT/${JC}_exres.log" 2>&1
  grep -q "ex res done" "$OUT/${JC}_exres.log" || { echo "$JC EXRES_FAILED"; exit 1; }
fi

# ---- pattern-recognition quality (BEFORE the reco file is deleted) ------------------------------
# The one check that decides whether a pitch difference may be read as a resolution difference:
# HDBSCAN's epsilon = 10 mm spans one pad at 8 x 12 mm and five at 2 mm, and a 2 mm track carries
# several times as many hits. If MERGED/SPLIT move between two cells of the matrix, the finder
# changed and the Ex comparison is not clean.
if ! grep -q "labelled hits placed in a cluster" "$OUT/${JC}_clu.log" 2>/dev/null; then
  echo "[$(date +%H:%M:%S)] $JC cluster eval"
  root -b -q -l "$SIM/cluster_eval_C14.C(\"$OUT/${JC}_reco.root\",\"$JC\",-1,20,0.10)" > "$OUT/${JC}_clu.log" 2>&1 || true
  mv -f "$SIM/diagnostics/cluster_${JC}.png" "$OUT/" 2>/dev/null || true
  mv -f "$SIM/cluster_${JC}.png" "$OUT/" 2>/dev/null || true
fi

rm -f "$OUT/${JC}_reco.root"
echo COMPLETED > "$OUT/$JC.marker"
echo "[$(date +%H:%M:%S)] $JC done: $(grep 'overall acceptance' "$OUT/${JC}_acc.log")"
