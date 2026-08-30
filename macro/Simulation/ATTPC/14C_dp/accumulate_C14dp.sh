#!/usr/bin/env bash
# ONE SAMPLE of the 14C(d,p)15C field x pad-pitch study: generate -> reco -> genfit -> acceptance
# + Ex resolution, for a single (level, field, pitch, seed).
#
#   ./accumulate_C14dp.sh <state> <B_T> <pad_mm> <seed> [nEvents]
#   ./accumulate_C14dp.sh gs 7.0 2.0 8021
#
# <state> is gs (the 1/2+ ground state) or ex0740 (the 5/2+ at 0.740 MeV). Those are the only two
# BOUND states of 15C -- S_n = 1.218 MeV -- so the pair is the whole bound spectrum, and 740 keV is
# a fair test of a resolution that runs 40-180 keV.
#
# WHY THIS CHANNEL. The (p,p') matrix found essentially no gain from pads or field, because its
# recoil protons are slow and come out near theta_lab 77 deg. A (d,p) proton at the angles where
# transfer peaks comes out BACKWARD (theta_lab 95-125 deg) with 3-8 MeV, where dEx/dKE is 1.2-2.1
# instead of 0.54 and dEx/dE_beam is ~0.047 instead of 0.004. The (p,p') campaign could not test
# that region at all: its protons never pass 90 deg.
#
# Everything except the target gas is held at the (p,p') campaign's configuration, so the two
# channels can be compared directly: same beam energy, same drift, same beam hole, same PSA and
# HDBSCAN settings, genfit with material effects + native CATIMA, chi2/ndf < 5 on GetKinematicsXtr.
#
# Resumable per stage, on evidence a stage finished rather than on a file existing.
set -eo pipefail

STATE=${1:?need a state: gs | ex0740}
BT=${2:?need the field in tesla}
PAD=${3:?need the pad pitch in mm, or -1 for the real AT-TPC pad plane}
SEED=${4:?need a seed}
NEV=${5:-8000}

# 15C has only TWO BOUND states -- the 1/2+ ground state and the 5/2+ at 0.740 -- because
# S_n = 1.218 MeV. The two higher entries are known but UNBOUND levels: a real 15C would emit a
# neutron. They are simulated anyway because the two-body proton kinematics are well defined and
# the proton is the only thing measured, so they are valid for a detector-response study. They are
# NOT a prediction of a bound spectrum.
case "$STATE" in
   gs)     EX=0.0   ;;   # 1/2+  bound
   ex0740) EX=0.740 ;;   # 5/2+  bound
   ex3103) EX=3.103 ;;   # unbound, detector response only
   ex4657) EX=4.657 ;;   # unbound, detector response only
   *) echo "unknown state '$STATE'"; exit 2 ;;
esac

REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
DP="$REPO/macro/Simulation/ATTPC/14C_dp"
PP="$REPO/macro/Simulation/ATTPC/14C_pp"      # run_reco_C14.C and acceptance_C14.C live here
HF="$PP/highfield"
UKF="$REPO/macro/Unpack_HDF5/a1954/UKF"
ROOTDIR=${DP_ROOT:-/mnt/f/a1954_C14dp_hf}

# the reaction, in amu -- the analysis macros take it as arguments now
MTGT=2.0141018      # d
MEJ=1.007825        # p
MRES=15.0105993     # 15C ground state
GEO=ATTPC_D300torr_v2

BTAG=$(awk -v b="$BT" 'BEGIN{printf "b%03d", b*100}')
PTAG=$(awk -v p="$PAD" 'BEGIN{ if (p>0) printf "%gmm", p; else printf "attpc" }')
# default measurement sigma for THIS pad plane (see the note at the genfit call)
MSDEF=$(awk -v p="$PAD" 'BEGIN{ print (p>0 && p<4) ? 0.35 : 0.6 }')
CFG="${BTAG}_${PTAG}"
SIMDIR="$ROOTDIR/sims_$BTAG"
OUT="$ROOTDIR/$CFG"
mkdir -p "$SIMDIR" "$OUT"

PAR="ATTPC.a1954_C14dp_${BTAG}.par"
[ -f "$REPO/parameters/$PAR" ] || { echo "MISSING PAR $PAR -- run make_dp_pars.sh"; exit 2; }

BKG=$(awk -v b="$BT" 'BEGIN{printf "%.1f", -10*b}')
BNEG=$(awk -v b="$BT" 'BEGIN{printf "%.2f", -b}')
J="${STATE}_s${SEED}"
JC="${J}_${CFG}"

set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"

echo "[cfg] $JC : Ex = $EX MeV, B = $BT T ($BKG kG), pads = $PAD mm, par = $PAR, gas $GEO, $NEV events"
[ -f "$OUT/$JC.marker" ] && { echo "$JC already COMPLETED"; exit 0; }

# ---- generation (shared per field) -------------------------------------------------------------
nsim=$(root -b -q -l -e "TFile*f=TFile::Open(\"$SIMDIR/${J}_sim.root\");TTree*t=(f&&!f->IsZombie())?(TTree*)f->Get(\"cbmsim\"):nullptr;printf(\"N %lld\\n\",t?t->GetEntries():-1);" 2>/dev/null | awk '/^N /{print $2}' || true)
nsim=${nsim:--1}
if [ "$nsim" -ne "$NEV" ]; then
  LOCK="$SIMDIR/${J}.lock"
  if mkdir "$LOCK" 2>/dev/null; then
    trap 'rmdir "$LOCK" 2>/dev/null || true' EXIT
    WORK="$SIMDIR/work_$J"; mkdir -p "$WORK/data"
    echo "[$(date +%H:%M:%S)] $J generating (had $nsim)"
    # CM_LO/CM_HI restrict the generated theta_cm range. The default is everything; setting them
    # concentrates the statistics where they are wanted, which for this channel is the BACKWARD
    # region (small theta_cm). The acceptance is a per-bin ratio and the resolution is measured per
    # slice, so neither is affected by the range -- only the density of events inside it is.
    ( cd "$WORK" && root -b -q -l "$DP/C14_dp_sim.C($NEV,${CM_LO:-2.0},${CM_HI:-178.0},\"TGeant4\",$BKG,\"$SIMDIR/${J}_sim.root\",$EX,$SEED,\"$GEO.root\")" ) > "$SIMDIR/${J}_gen.log" 2>&1
    # Compare the range NUMERICALLY. ROOT prints "2 - 60" for arguments 2.0 and 60.0, so a string
    # match on the argument as written fails on a run that is perfectly correct. This is the exact
    # trap the 46Ar campaign hit on its excitation energy, written down at the time, and repeated
    # here anyway -- which is the argument for never checking a numeric argument as a string.
    cmlog=$(awk '/CM angular range:/{print $4, $6; exit}' "$SIMDIR/${J}_gen.log")
    awk -v got="$cmlog" -v a="${CM_LO:-2.0}" -v b="${CM_HI:-178.0}" \
        'BEGIN{n=split(got,g," "); ok=(n==2 && g[1]==g[1]+0 && (g[1]-a<1e-6 && a-g[1]<1e-6) && (g[2]-b<1e-6 && b-g[2]<1e-6)); exit !ok}' \
        || { echo "$J CM_RANGE_NOT_APPLIED (log says '"'"'${cmlog:-none}'"'"', wanted ${CM_LO:-2.0} - ${CM_HI:-178.0})"; exit 1; }
    grep -q "$GEO" "$SIMDIR/${J}_gen.log" || { echo "$J WRONG_GAS"; exit 1; }
    grep -q "RNG seed requested: $SEED" "$SIMDIR/${J}_gen.log" || { echo "$J SEED_NOT_APPLIED"; exit 1; }
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

cd "$PP"

# ---- digitisation + reconstruction --------------------------------------------------------------
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

# ---- genfit --------------------------------------------------------------------------------------
if [ -s "$OUT/${JC}_genfit.root" ] && grep -q "dE/dx from CATIMA" "$OUT/${JC}_fit.log" 2>/dev/null; then
  echo "[$(date +%H:%M:%S)] $JC fit already complete"
else
  echo "[$(date +%H:%M:%S)] $JC genfit at B = $BNEG T, gas $GEO"
  # backwardSeedFix = kTRUE (argument 14), and it is NOT optional in this channel. AtGenfitter
  # seeds a track from its lowest-z end with the momentum pointing to +z, which REFLECTS a
  # backward-going track into the forward hemisphere. The (p,p') campaign never noticed because
  # its recoil protons never pass theta_lab 90 deg; here the whole point is the ones that do. The
  # first (d,p) run left it kFALSE and got ZERO reconstructed protons below theta_cm 45 deg.
  # MEASSIGMA IS PER PAD PLANE, not global. chi2/ndf implies the real hit residuals are
  # 0.59-0.64 mm on the AT-TPC plane and 0.32-0.35 mm on the 2 mm plane -- finer pads really do
  # localise better, so one value cannot be right for both. Measured from
  # measSigma*sqrt(chi2/ndf) per configuration (dp_perbin_C14.C, panel B).
  # Previously (0.6 mm everywhere): 0.6 mm, not the 4.0 inherited from (p,p'). Measured -- chi2/ndf
  # scales as 1/measSigma^2 (0.015 / 0.055 / 0.207 / 0.548 / 2.111 at 4.0 / 2.0 / 1.0 / 0.6 / 0.3),
  # so the real fit residuals are ~0.45 mm and 4 mm left chi2/ndf at 0.015, making the quality cut
  # meaningless. 0.6 puts chi2/ndf at 0.55 and gives the best energy spread (1.00 % -> 0.84 %).
  root -b -q -l "$UKF/pipeline/fitGenfit_C14.C(\"$JC\",-1,\"$OUT/\",\"\",\"$OUT/\",$BNEG,2,5,\"\",${MEASSIGMA:-$MSDEF},10.0,170.0,kTRUE,kTRUE,\"proton\",\"$GEO\",kFALSE,kTRUE,0.0,1,kTRUE,kFALSE,kFALSE,kFALSE)" > "$OUT/${JC}_fit.log" 2>&1
  [ -s "$OUT/${JC}_genfit.root" ] || { echo "$JC FIT_FAILED"; exit 1; }
  grep -q "dE/dx from CATIMA" "$OUT/${JC}_fit.log" || { echo "$JC CATIMA_NOT_ENABLED"; exit 1; }
fi

# ---- acceptance -----------------------------------------------------------------------------------
if ! grep -q "overall acceptance" "$OUT/${JC}_acc.log" 2>/dev/null; then
  echo "[$(date +%H:%M:%S)] $JC acceptance"
  root -b -q -l "$PP/acceptance_C14.C(\"$SIMDIR/${J}_sim.root\",\"$OUT/${JC}_genfit.root\",\"$JC\",$EX,159.75,5.0,36,180.0,10.0,0.5,2.0,kTRUE,-1e9,1e9,$MTGT,$MEJ,$MRES)" > "$OUT/${JC}_acc.log" 2>&1
  grep -q "overall acceptance" "$OUT/${JC}_acc.log" || { echo "$JC ACC_FAILED"; exit 1; }
  mv -f "$PP/diagnostics/acceptance_${JC}.root" "$OUT/" 2>/dev/null || true
  rm -f "$PP/diagnostics/acceptance_${JC}.png"
fi

# ---- excitation-energy resolution -------------------------------------------------------------------
if ! grep -q "ex res done" "$OUT/${JC}_exres.log" 2>/dev/null; then
  echo "[$(date +%H:%M:%S)] $JC Ex resolution"
  root -b -q -l "$HF/ex_res_C14_hf.C(\"$SIMDIR/${J}_sim.root\",\"$OUT/${JC}_genfit.root\",\"$JC\",$EX,159.75,5.0,kTRUE,\"$OUT/\",10.0,0.5,2.0,$MTGT,$MEJ,$MRES)" > "$OUT/${JC}_exres.log" 2>&1
  grep -q "ex res done" "$OUT/${JC}_exres.log" || { echo "$JC EXRES_FAILED"; exit 1; }
fi

# ---- pattern-recognition quality, before the reco file goes ------------------------------------------
if ! grep -q "labelled hits placed in a cluster" "$OUT/${JC}_clu.log" 2>/dev/null; then
  root -b -q -l "$PP/cluster_eval_C14.C(\"$OUT/${JC}_reco.root\",\"$JC\",-1,20,0.10)" > "$OUT/${JC}_clu.log" 2>&1 || true
  mv -f "$PP/diagnostics/cluster_${JC}.png" "$OUT/" 2>/dev/null || true
fi

rm -f "$OUT/${JC}_reco.root"
echo COMPLETED > "$OUT/$JC.marker"
echo "[$(date +%H:%M:%S)] $JC done: $(grep 'overall acceptance' "$OUT/${JC}_acc.log")"
