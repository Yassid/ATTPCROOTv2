#!/usr/bin/env bash
# ONE SAMPLE of the 10Be(t,p)12Be field x pad-pitch study: generate -> reco -> genfit -> acceptance
# + Ex resolution, for a single (level, field, pitch, seed).
#
#   ./accumulate_Be10tp.sh <state> <B_T> <pad_mm> <seed> [nEvents]
#   ./accumulate_Be10tp.sh gs 7.0 2.0 9021
#
# <state> is one of the four BOUND levels of 12Be (S_n = 3.171 MeV, so all four are bound --
# unlike 15C, where the (d,p) study had only two):
#
#   gs       0.000  0+
#   ex2109   2.109  2+
#   ex2251   2.251  0+_2, the intruder      <- populated 5x more weakly than the others
#   ex2715   2.715  1-
#
# THE 5x SUPPRESSION OF THE 0+_2 IS APPLIED AT COMBINATION TIME, NOT HERE. Every level is
# simulated at FULL statistics so that its acceptance and its Ex resolution are each measured to
# the same precision; the relative population enters only when the four are summed into one
# spectrum (tp_spectrum_Be10.C, weight 0.2 on ex2251). Generating a fifth of the events instead
# would have made the level's own acceptance the worst-measured of the four, which is backwards --
# it is the one whose extraction is hardest.
#
# Everything except the target gas and the reaction is held at the 14C(d,p) campaign's
# configuration, so the two channels can be compared directly: same beam energy PER NUCLEON
# (11.5 MeV/u), same 300 torr, same drift, same beam hole, same PSA and HDBSCAN settings, genfit
# with material effects + native CATIMA, chi2/ndf < 5 on GetKinematicsXtr.
#
# Resumable per stage, on evidence a stage finished rather than on a file existing.
set -eo pipefail

STATE=${1:?need a state: gs | ex2109 | ex2251 | ex2715}
BT=${2:?need the field in tesla}
PAD=${3:?need the pad pitch in mm, or -1 for the real AT-TPC pad plane}
SEED=${4:?need a seed}
NEV=${5:-16000}

case "$STATE" in
   gs)     EX=0.0   ;;   # 0+   ground state
   ex2109) EX=2.109 ;;   # 2+
   ex2251) EX=2.251 ;;   # 0+_2 intruder, only 142 keV above the 2+
   ex2715) EX=2.715 ;;   # 1-
   *) echo "unknown state '$STATE'"; exit 2 ;;
esac

REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
TP="$REPO/macro/Simulation/ATTPC/10Be_tp"
PP="$REPO/macro/Simulation/ATTPC/14C_pp"      # run_reco_C14.C and acceptance_C14.C live here
HF="$PP/highfield"
UKF="$REPO/macro/Unpack_HDF5/a1954/UKF"
ROOTDIR=${TP_ROOT:-/mnt/f/Be10_tp}

# the reaction, in amu -- the analysis macros take every mass as an argument
MBEAM=10.0135341    # 10Be
MTGT=3.0160493      # t
MEJ=1.007825        # p
MRES=12.0269221     # 12Be ground state
GEO=ATTPC_T300torr
# Beam energy AT THE REACTION VERTEX, measured from MC truth by check_vertex_beam_Be10.C:
# 115.0 MeV at the window, mean vertex at z = 503 mm, 2.80 MeV lost getting there -> 112.20.
# This is the constant the two-body inversion in acceptance/ex_res must use; using 115 would put
# a coherent offset on every reconstructed Ex.
EBEAM=112.20

BTAG=$(awk -v b="$BT" 'BEGIN{printf "b%03d", b*100}')
PTAG=$(awk -v p="$PAD" 'BEGIN{ if (p>0) printf "%gmm", p; else printf "attpc" }')
# default measurement sigma for THIS pad plane. Measured on the (d,p) campaign from
# measSigma*sqrt(chi2/ndf): the real hit residuals are 0.59-0.64 mm on the AT-TPC plane and
# 0.32-0.35 mm on a 2 mm plane. One global value is wrong for half the matrix.
MSDEF=$(awk -v p="$PAD" 'BEGIN{ print (p>0 && p<4) ? 0.35 : 0.6 }')
CFG="${BTAG}_${PTAG}"
SIMDIR="$ROOTDIR/sims_$BTAG"
OUT="$ROOTDIR/$CFG"
mkdir -p "$SIMDIR" "$OUT"

PAR="ATTPC.Be10tp_${BTAG}.par"
[ -f "$REPO/parameters/$PAR" ] || { echo "MISSING PAR $PAR -- run make_tp_pars.sh"; exit 2; }

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
    ( cd "$WORK" && root -b -q -l "$TP/Be10_tp_sim.C($NEV,${CM_LO:-2.0},${CM_HI:-178.0},\"TGeant4\",$BKG,\"$SIMDIR/${J}_sim.root\",$EX,$SEED,\"$GEO.root\")" ) > "$SIMDIR/${J}_gen.log" 2>&1
    # Compare the CM range NUMERICALLY: ROOT prints "2 - 60" for the arguments 2.0 and 60.0, so a
    # string match on the argument as written fails on a run that is perfectly correct.
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
  # backwardSeedFix = kTRUE (argument 14), and it is NOT optional here. AtGenfitter seeds a track
  # from its lowest-z end with the momentum pointing to +z, which REFLECTS a backward-going track
  # into the forward hemisphere. The (d,p) campaign's first run left it kFALSE and got ZERO
  # reconstructed protons below theta_cm 45 deg -- which read as physics.
  root -b -q -l "$UKF/pipeline/fitGenfit_C14.C(\"$JC\",-1,\"$OUT/\",\"\",\"$OUT/\",$BNEG,2,5,\"\",${MEASSIGMA:-$MSDEF},10.0,170.0,kTRUE,kTRUE,\"proton\",\"$GEO\",kFALSE,kTRUE,0.0,1,kTRUE,kFALSE,kFALSE,kFALSE)" > "$OUT/${JC}_fit.log" 2>&1
  [ -s "$OUT/${JC}_genfit.root" ] || { echo "$JC FIT_FAILED"; exit 1; }
  grep -q "dE/dx from CATIMA" "$OUT/${JC}_fit.log" || { echo "$JC CATIMA_NOT_ENABLED"; exit 1; }
fi

# ---- acceptance -----------------------------------------------------------------------------------
if ! grep -q "overall acceptance" "$OUT/${JC}_acc.log" 2>/dev/null; then
  echo "[$(date +%H:%M:%S)] $JC acceptance"
  root -b -q -l "$PP/acceptance_C14.C(\"$SIMDIR/${J}_sim.root\",\"$OUT/${JC}_genfit.root\",\"$JC\",$EX,$EBEAM,5.0,36,180.0,10.0,0.5,2.0,kTRUE,-1e9,1e9,$MTGT,$MEJ,$MRES,$MBEAM)" > "$OUT/${JC}_acc.log" 2>&1
  grep -q "overall acceptance" "$OUT/${JC}_acc.log" || { echo "$JC ACC_FAILED"; exit 1; }
  mv -f "$PP/diagnostics/acceptance_${JC}.root" "$OUT/" 2>/dev/null || true
  rm -f "$PP/diagnostics/acceptance_${JC}.png"
fi

# ---- excitation-energy resolution -------------------------------------------------------------------
if ! grep -q "ex res done" "$OUT/${JC}_exres.log" 2>/dev/null; then
  echo "[$(date +%H:%M:%S)] $JC Ex resolution"
  root -b -q -l "$HF/ex_res_C14_hf.C(\"$SIMDIR/${J}_sim.root\",\"$OUT/${JC}_genfit.root\",\"$JC\",$EX,$EBEAM,5.0,kTRUE,\"$OUT/\",10.0,0.5,2.0,$MTGT,$MEJ,$MRES,$MBEAM)" > "$OUT/${JC}_exres.log" 2>&1
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
