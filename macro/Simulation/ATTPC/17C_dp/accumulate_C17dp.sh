#!/usr/bin/env bash
# ONE SAMPLE of the 17C(d,p)18C proposal simulation: generate -> reco -> genfit -> acceptance
# + Ex resolution, for a single 18C level.
#
#   ./accumulate_C17dp.sh <state> <seed> [nEvents]
#   ./accumulate_C17dp.sh gs 9001 12000
#
# CHANNEL. 17C at 8.37 MeV/u = 142.29 MeV on D2 at 300 torr, B = 2.85 T, real AT-TPC pad plane --
# the SOLARIS + AT-TPC configuration of the 17C M_n/M_p proposal, during its deuterium day.
# Q(17C(d,p)18C) = +1.959 MeV.
#
# <state> selects the 18C level. 18C has S_n = 4184 keV and FOUR bound states, so this list is the
# whole bound spectrum -- there is no unbound-level padding as there was in 14C(d,p):
#     gs      0        0+          ground state
#     ex1588  1.588    2+          15.5 ps, B(E2) = 0.000364 e2b2; the state the proposal discusses
#     ex2515  2.515    (2+)        < 3.2 ps
#     ex3972  3.972    (2,3)+
# The 927 keV gap between the 1588 and the 2515 is the separation that has to be resolved, against
# the ~300 keV the proposal quotes for the AT-TPC.
#
# WHY THE SETTINGS ARE NOT TUNED HERE. Everything except the beam and the reaction is held at the
# 14C(d,p) reference campaign's configuration (macro/Simulation/ATTPC/14C_dp/RESULTS.md), which was
# debugged against data: same gas, same field, same drift, same beam hole, same PSA and HDBSCAN
# settings, genfit with material effects + native CATIMA, chi2/ndf < 5 on GetKinematicsXtr. A
# difference between the two channels is then the reaction, not a setting.
#
# Resumable per stage, on evidence a stage finished rather than on a file existing.
set -eo pipefail

STATE=${1:?need a state: gs | ex1588 | ex2515 | ex3972}
SEED=${2:?need a seed}
NEV=${3:-12000}

case "$STATE" in
   gs)     EX=0.0   ;;   # 0+      bound
   ex1588) EX=1.588 ;;   # 2+      bound
   ex2515) EX=2.515 ;;   # (2+)    bound
   ex3972) EX=3.972 ;;   # (2,3)+  bound
   *) echo "unknown state '$STATE' (gs | ex1588 | ex2515 | ex3972)"; exit 2 ;;
esac

REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
DP="$REPO/macro/Simulation/ATTPC/17C_dp"
PP="$REPO/macro/Simulation/ATTPC/14C_pp"      # run_reco_C14.C and acceptance_C14.C live here
HF="$PP/highfield"                            # ex_res_C14_hf.C
UKF="$REPO/macro/Unpack_HDF5/a1954/UKF"       # fitGenfit_C14.C
ROOTDIR=${DP_ROOT:-/mnt/f/C17dp}

# THE ANALYSIS MACROS ARE BEAM-AGNOSTIC -- acceptance_C14.C and ex_res_C14_hf.C take the target,
# ejectile, residual AND beam masses as trailing arguments (the last was added 2026-08-31 for
# 10Be(t,p)). Nothing about "C14" in their names is baked into the physics; only the defaults are.
MTGT=2.0141018      # d
MEJ=1.007825        # p
MRES=18.0267519     # 18C ground state, AME2020
MBEAM=17.0225787    # 17C, AME2020
GEO=ATTPC_D300torr_v2

# B = 2.85 T, in the DATA (negative) convention -- see the sign note in C17_dp_sim.C.
BT=2.85
BKG=-28.5           # kG, for the transport field
BNEG=-2.85          # T, for genfit
PAD=-1              # the real AT-TPC pad plane
MEASSIGMA=${MEASSIGMA:-0.6}   # mm. MEASURED for this pad plane in 14C(d,p): chi2/ndf implies real
                              # hit residuals of 0.59-0.64 mm on the AT-TPC plane. 4.0 (inherited
                              # from the (p,p') campaign) leaves chi2/ndf at 0.015 and makes the
                              # quality cut meaningless.

# BEAM ENERGY for the constant-Ebeam Ex reconstruction. 142.29 MeV entering, 127.49 leaving the
# metre (CATIMA in this gas), and the vertex is uniform in z, so the mean beam energy AT THE VERTEX
# is 135.0 MeV. That is what a constant-Ebeam analysis should use; a vertex-dependent correction
# does better still (in 14C(d,p) it took backward sigma(Ex) from 0.178 to 0.064 MeV) and is a
# separate, software-only step.
EBEAM=${EBEAM:-135.0}

CFG=b285_attpc
SIMDIR="$ROOTDIR/sims"
OUT="$ROOTDIR/$CFG"
mkdir -p "$SIMDIR" "$OUT"

PAR="ATTPC.C17dp_D300torr_b285.par"
[ -f "$REPO/parameters/$PAR" ] || { echo "MISSING PAR $PAR -- run make_c17dp_par.sh"; exit 2; }

J="${STATE}_s${SEED}"
JC="${J}_${CFG}"

set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"

echo "[cfg] $JC : Ex = $EX MeV, B = $BT T ($BKG kG), AT-TPC pads, par = $PAR, gas $GEO, $NEV events"
[ -f "$OUT/$JC.marker" ] && { echo "$JC already COMPLETED"; exit 0; }

# Entries in a ROOT file's cbmsim tree, or -1 if the file is missing, unopenable, or was never
# closed. THIS IS THE ONLY ACCEPTABLE COMPLETION TEST FOR A STAGE THAT WRITES A ROOT FILE.
#
# The campaign's first run failed here and it is worth the comment. The genfit resume check used to
# be "the file is non-empty AND the log says CATIMA", which is satisfied the moment genfit STARTS
# writing. A second campaign driver was accidentally left running on this same DP_ROOT, so while
# driver A was still writing ex2515_genfit.root, driver B read "fit already complete" and handed
# the half-written file to the acceptance stage, which got "probably not closed / missing cbmsim"
# and failed. A check keyed on a file EXISTING, or on a log line that is printed before the work is
# done, is not evidence the work finished -- only reading the finished product back is.
#
# `|| true` IS LOAD-BEARING, do not tidy it away. This script runs under `set -eo pipefail`, and
# `root -b -q -l -e ...` exits NON-ZERO (8) even when it printed the answer perfectly. Without the
# guard, pipefail propagates that through the pipe, the command substitution fails, and `set -e`
# kills the whole sample silently -- which is exactly what happened on the second campaign run: all
# four samples died immediately after the generation check with no error message and no genfit log,
# having just spent two hours on reco. The generation check that this helper was factored out of
# always had the guard; the extracted copy lost it.
nEntries() {
   local f="$1" n
   [ -s "$f" ] || { echo -1; return 0; }
   n=$(root -b -q -l -e "TFile*f=TFile::Open(\"$f\");TTree*t=(f&&!f->IsZombie())?(TTree*)f->Get(\"cbmsim\"):nullptr;printf(\"N %lld\\n\",t?t->GetEntries():-1);" 2>/dev/null \
      | awk '/^N /{print $2; exit}' || true)
   echo "${n:--1}"
   return 0
}

# ---- generation --------------------------------------------------------------------------------
nsim=$(root -b -q -l -e "TFile*f=TFile::Open(\"$SIMDIR/${J}_sim.root\");TTree*t=(f&&!f->IsZombie())?(TTree*)f->Get(\"cbmsim\"):nullptr;printf(\"N %lld\\n\",t?t->GetEntries():-1);" 2>/dev/null | awk '/^N /{print $2}' || true)
nsim=${nsim:--1}
if [ "$nsim" -ne "$NEV" ]; then
  LOCK="$SIMDIR/${J}.lock"
  if mkdir "$LOCK" 2>/dev/null; then
    trap 'rmdir "$LOCK" 2>/dev/null || true' EXIT
    WORK="$SIMDIR/work_$J"; mkdir -p "$WORK/data"
    echo "[$(date +%H:%M:%S)] $J generating (had $nsim)"
    # CM_LO/CM_HI restrict the generated theta_cm range; the default is everything. The generator
    # is flat in cos(theta_cm), and acceptance is a per-bin ratio while resolution is measured per
    # slice, so narrowing the range changes only the density of events inside it.
    ( cd "$WORK" && root -b -q -l "$DP/C17_dp_sim.C($NEV,${CM_LO:-2.0},${CM_HI:-178.0},\"TGeant4\",$BKG,\"$SIMDIR/${J}_sim.root\",$EX,$SEED,\"$GEO.root\")" ) > "$SIMDIR/${J}_gen.log" 2>&1
    # Compare the CM range NUMERICALLY. ROOT prints "2 - 178" for arguments 2.0 and 178.0, so a
    # string match on the argument as written fails on a run that is perfectly correct. This is the
    # trap the 46Ar campaign hit on its excitation energy and 14C(d,p) hit again on this same line.
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
nreco=$(nEntries "$OUT/${JC}_reco.root")
if [ "$nreco" -eq "$NEV" ] && grep -q "sim reco done" "$OUT/${JC}_reco.log" 2>/dev/null; then
  echo "[$(date +%H:%M:%S)] $JC reco already complete ($nreco entries)"
else
  echo "[$(date +%H:%M:%S)] $JC reco (had $nreco)"
  root -b -q -l "run_reco_C14.C(\"$SIMDIR/${J}_sim.root\",\"$OUT/${JC}_reco.root\",\"$PAR\",20,20,8,0,30.0,\"mover\",$PAD,500.0,kFALSE)" > "$OUT/${JC}_reco.log" 2>&1
  grep -q "sim reco done" "$OUT/${JC}_reco.log" || { echo "$JC RECO_FAILED"; exit 1; }
  nreco=$(nEntries "$OUT/${JC}_reco.root")
  [ "$nreco" -eq "$NEV" ] || { echo "$JC RECO_INCOMPLETE ($nreco of $NEV)"; exit 1; }
fi

# ---- genfit --------------------------------------------------------------------------------------
ngf=$(nEntries "$OUT/${JC}_genfit.root")
# The fit output must READ BACK with the same number of entries as its input. acceptance_C14.C and
# ex_res_C14_hf.C both index the sim and fit trees by entry number and refuse to run if the counts
# disagree, so anything short here is not merely incomplete, it is unusable.
if [ "$ngf" -eq "$NEV" ] && grep -q "dE/dx from CATIMA" "$OUT/${JC}_fit.log" 2>/dev/null; then
  echo "[$(date +%H:%M:%S)] $JC fit already complete ($ngf entries)"
else
  echo "[$(date +%H:%M:%S)] $JC genfit at B = $BNEG T, gas $GEO"
  # backwardSeedFix = kTRUE (argument 14) is NOT optional in this channel. AtGenfitter seeds a
  # track from its lowest-z end with the momentum pointing to +z, which REFLECTS a backward-going
  # track into the forward hemisphere. The transfer peak here is backward in the lab, exactly where
  # that bites: the first 14C(d,p) run left it kFALSE and got ZERO protons below theta_cm 45 deg.
  root -b -q -l "$UKF/pipeline/fitGenfit_C14.C(\"$JC\",-1,\"$OUT/\",\"\",\"$OUT/\",$BNEG,2,5,\"\",$MEASSIGMA,10.0,170.0,kTRUE,kTRUE,\"proton\",\"$GEO\",kFALSE,kTRUE,0.0,1,kTRUE,kFALSE,kFALSE,kFALSE)" > "$OUT/${JC}_fit.log" 2>&1
  grep -q "dE/dx from CATIMA" "$OUT/${JC}_fit.log" || { echo "$JC CATIMA_NOT_ENABLED"; exit 1; }
  ngf=$(nEntries "$OUT/${JC}_genfit.root")
  [ "$ngf" -eq "$NEV" ] || { echo "$JC FIT_INCOMPLETE ($ngf of $NEV) -- the file did not close"; exit 1; }
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
