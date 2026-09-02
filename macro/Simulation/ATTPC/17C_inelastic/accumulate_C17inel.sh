#!/usr/bin/env bash
# ONE SAMPLE of the 17C M_n/M_p proposal simulation: generate -> reco -> genfit -> acceptance
# + Ex resolution, for a single (channel, level, field).
#
#   ./accumulate_C17inel.sh <channel> <state> <field_T> <seed> [nEvents]
#   ./accumulate_C17inel.sh pp gs   2.85 5101 8000
#   ./accumulate_C17inel.sh dd ex332 4.0  5212 8000
#
# CHANNEL. 17C at 8.37 MeV/u = 142.29 MeV (ReA6, 940 pps) on 300 torr H2 for (p,p') or 300 torr D2
# for (d,d') -- the two days of the proposal C17p_FRIB_Proposal.pdf. The AT-TPC is both target and
# tracker; the fitted track is the light recoil.
#
# <state> selects the 17C level. S_n(17C) = 733 keV, so these three ARE the bound spectrum:
#     gs     0        3/2+
#     ex217  0.217    1/2+     proposal target
#     ex332  0.332    5/2+     proposal target   (ENSDF adopts 331 keV)
# The two excited states are 115 keV apart and the ground state is 217 keV below the 1/2+, so all
# three live inside one 300 keV resolution width with the elastic far stronger than either
# inelastic. That is the measurement, and quantifying it is what this campaign is for.
#
# <field_T> is 2.85 (the field Ref.[24] ran at, and this campaign's nominal) or 4.0 (SOLARIS's
# design field). The kinematics gate -- inel_kinematics_C17.C -- says the field is the ONLY lever
# on this doublet: sigma(Ex) is flat in angle and set entirely by sigma(KE), which the 14C matrix
# measured as 0.343 MeV at 2.85 T and 0.067 MeV at 4 T on this same pad plane. Hence a field axis
# and no pad-pitch axis: this proposal runs the conventional AT-TPC pad plane.
#
# WHY THE SETTINGS ARE NOT TUNED HERE. Everything except the beam, the target and the field is held
# at the 14C(p,p') reference campaign's configuration (14C_pp/highfield/), which was debugged
# against a1954 data: same drift, same beam hole, same PSA and HDBSCAN settings, genfit with
# material effects + native CATIMA, chi2/ndf < 5 on GetKinematicsXtr. A difference against that
# campaign is then the beam and the reaction, not a setting.
#
# Resumable per stage, on evidence a stage finished rather than on a file existing.
set -eo pipefail

CHAN=${1:?need a channel: pp | dd}
STATE=${2:?need a state: gs | ex217 | ex332}
BT=${3:?need a field in T: 2.85 | 4.0}
SEED=${4:?need a seed}
NEV=${5:-8000}

case "$STATE" in
   gs)    EX=0.0   ;;   # 3/2+  ground state, elastic
   ex217) EX=0.217 ;;   # 1/2+
   ex332) EX=0.332 ;;   # 5/2+
   *) echo "unknown state '$STATE' (gs | ex217 | ex332)"; exit 2 ;;
esac

REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
IN="$REPO/macro/Simulation/ATTPC/17C_inelastic"
PP="$REPO/macro/Simulation/ATTPC/14C_pp"      # run_reco_C14.C, acceptance_C14.C, cluster_eval_C14.C
HF="$PP/highfield"                            # ex_res_C14_hf.C
UKF="$REPO/macro/Unpack_HDF5/a1954/UKF"       # fitGenfit_C14.C
ROOTDIR=${INEL_ROOT:-/media/yassid/Seagate Hub/ATTPC/C17_inel}

# THE ANALYSIS MACROS ARE BEAM-AGNOSTIC -- acceptance_C14.C and ex_res_C14_hf.C take the target,
# ejectile, residual AND beam masses as trailing arguments. Nothing about "C14" in their names is
# baked into the physics; only the defaults are.
# For inelastic scattering the residual IS the beam: 17C, excited by EX.
MBEAM=17.0225787    # 17C, AME2020
MRES=17.0225787     # 17C again -- this is scattering, not transfer
case "$CHAN" in
   pp) MTGT=1.007825;   MEJ=1.007825;   GEO=ATTPC_H300torr_RT; PART=proton;   EJPDG=2212 ;;
   dd) MTGT=2.0141018;  MEJ=2.0141018;  GEO=ATTPC_D300torr_v2; PART=deuteron; EJPDG=1000010020 ;;
   *) echo "unknown channel '$CHAN' (pp | dd)"; exit 2 ;;
esac

# Field, in the DATA (negative) convention -- see the sign note in C17_inel_sim.C. BTAG is the
# 3-digit centi-tesla code the 14C campaign uses, so par names and output dirs line up with it.
BTAG=$(awk -v b="$BT" 'BEGIN{printf "b%03d", b*100}')
BKG=$(awk -v b="$BT" 'BEGIN{printf "%.1f", -b*10}')    # kG, transport field
BNEG=$(awk -v b="$BT" 'BEGIN{printf "%.2f", -b}')      # T, genfit
PAD=-1              # the real AT-TPC pad plane -- this proposal has no pad-plane development

# Par files. All four already exist in parameters/ and are gas- and field-matched; CoefT at 4 T is
# the Magboltz-scaled value the 14C high-field campaign derived, anchored on the a1954-tuned 9e-4
# at 2.85 T. Nothing new is generated here, which is why there is no make_*_par.sh in this
# directory.
case "${CHAN}_${BTAG}" in
   pp_b285) PAR="ATTPC.a1954_C14_hf_b285.par"    ;;   # H2 300 torr, 2.85 T
   pp_b400) PAR="ATTPC.a1954_C14_hf_b400.par"    ;;   # H2 300 torr, 4.00 T
   dd_b285) PAR="ATTPC.C17dp_D300torr_b285.par"  ;;   # D2 300 torr, 2.85 T
   dd_b400) PAR="ATTPC.a1954_C14dp_b400.par"     ;;   # D2 300 torr, 4.00 T
   *) echo "no par for ${CHAN}_${BTAG}"; exit 2 ;;
esac
[ -f "$REPO/parameters/$PAR" ] || { echo "MISSING PAR $PAR"; exit 2; }

# genfit's assumed per-hit measurement error. 4.0 mm is what the 14C(p,p') reference campaign used,
# and it is kept here so this campaign's resolution numbers are comparable with the sigma(KE)
# values that the kinematics gate extrapolates from. It is NOT the measured hit residual: the
# (d,p) arm found 0.59-0.64 mm on this pad plane, at which the chi2/ndf < 5 cut actually bites
# rather than passing everything. Run one cell with MEASSIGMA=0.6 as a control before quoting a
# tail fraction.
MEASSIGMA=${MEASSIGMA:-4.0}

# BEAM ENERGY for the constant-Ebeam Ex reconstruction. 142.29 MeV entering, ~127.5 leaving the
# metre, vertex uniform in z -> mean beam energy AT THE VERTEX is 135.0 MeV. A vertex-dependent
# correction does better (on 14C(p,p') it took sigma(Ex) from 0.265 to 0.176 MeV) and is a
# separate, software-only step applied in inel_summary_C17.C.
EBEAM=${EBEAM:-135.0}

CFG="${CHAN}_${BTAG}"
SIMDIR="$ROOTDIR/sims"
OUT="$ROOTDIR/$CFG"
mkdir -p "$SIMDIR" "$OUT"

# The generation depends on the channel (gas), the field (it bends the tracks) and the level, so
# every one of those is in the key. Two cells never share a sim file.
J="${CHAN}_${STATE}_${BTAG}_s${SEED}"
JC="$J"

set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"

echo "[cfg] $JC : Ex = $EX MeV, B = $BT T ($BKG kG), AT-TPC pads, par = $PAR, gas $GEO, fit as $PART, $NEV events"
[ -f "$OUT/$JC.marker" ] && { echo "$JC already COMPLETED"; exit 0; }

# Entries in a ROOT file's cbmsim tree, or -1 if the file is missing, unopenable, or was never
# closed. THIS IS THE ONLY ACCEPTABLE COMPLETION TEST FOR A STAGE THAT WRITES A ROOT FILE: a check
# keyed on a file EXISTING, or on a log line printed before the work is done, is not evidence the
# work finished. The (d,p) arm lost two hours to exactly that.
#
# `|| true` IS LOAD-BEARING, do not tidy it away. Under `set -eo pipefail`, `root -b -q -l -e ...`
# exits NON-ZERO (8) even when it printed the answer perfectly; without the guard pipefail
# propagates that and `set -e` kills the sample silently.
nEntries() {
   local f="$1" n
   [ -s "$f" ] || { echo -1; return 0; }
   n=$(root -b -q -l -e "TFile*f=TFile::Open(\"$f\");TTree*t=(f&&!f->IsZombie())?(TTree*)f->Get(\"cbmsim\"):nullptr;printf(\"N %lld\\n\",t?t->GetEntries():-1);" 2>/dev/null \
      | awk '/^N /{print $2; exit}' || true)
   echo "${n:--1}"
   return 0
}

# ---- generation --------------------------------------------------------------------------------
nsim=$(nEntries "$SIMDIR/${J}_sim.root")
if [ "$nsim" -ne "$NEV" ]; then
  LOCK="$SIMDIR/${J}.lock"
  if mkdir "$LOCK" 2>/dev/null; then
    trap 'rmdir "$LOCK" 2>/dev/null || true' EXIT
    WORK="$SIMDIR/work_$J"; mkdir -p "$WORK/data"
    echo "[$(date +%H:%M:%S)] $J generating (had $nsim)"
    # CM_LO/CM_HI restrict the generated theta_cm range; the default starts at 10 deg because below
    # that the recoil carries under 220 keV and cannot make a track. The generator is flat in
    # cos(theta_cm), and acceptance is a per-bin ratio while resolution is measured per slice, so
    # narrowing the range changes only the density of events inside it.
    ( cd "$WORK" && root -b -q -l "$IN/C17_inel_sim.C($NEV,\"$CHAN\",$EX,$SEED,$BKG,\"$SIMDIR/${J}_sim.root\",${CM_LO:-10.0},${CM_HI:-178.0},\"TGeant4\")" ) > "$SIMDIR/${J}_gen.log" 2>&1
    # Compare the CM range NUMERICALLY. ROOT prints "10 - 178" for arguments 10.0 and 178.0, so a
    # string match on the argument as written fails on a run that is perfectly correct. This trap
    # has now bitten the 46Ar, 14C(d,p) and 17C(d,p) campaigns on this same line.
    cmlog=$(awk '/CM angular range:/{print $4, $6; exit}' "$SIMDIR/${J}_gen.log")
    awk -v got="$cmlog" -v a="${CM_LO:-10.0}" -v b="${CM_HI:-178.0}" \
        'BEGIN{n=split(got,g," "); ok=(n==2 && g[1]==g[1]+0 && (g[1]-a<1e-6 && a-g[1]<1e-6) && (g[2]-b<1e-6 && b-g[2]<1e-6)); exit !ok}' \
        || { echo "$J CM_RANGE_NOT_APPLIED (log says '"'"'${cmlog:-none}'"'"', wanted ${CM_LO:-10.0} - ${CM_HI:-178.0})"; exit 1; }
    grep -q "$GEO" "$SIMDIR/${J}_gen.log" || { echo "$J WRONG_GAS"; exit 1; }
    grep -q "channel $CHAN " "$SIMDIR/${J}_gen.log" || { echo "$J WRONG_CHANNEL"; exit 1; }
    grep -q "RNG seed requested: $SEED" "$SIMDIR/${J}_gen.log" || { echo "$J SEED_NOT_APPLIED"; exit 1; }
    rm -rf "$WORK"
    rmdir "$LOCK" 2>/dev/null || true; trap - EXIT
  else
    echo "[$(date +%H:%M:%S)] $J generation held by another job, waiting"
    for _ in $(seq 1 720); do sleep 10; [ -d "$LOCK" ] || break; done
    [ -d "$LOCK" ] && { echo "$J GEN_LOCK_TIMEOUT"; exit 1; }
  fi
  nsim=$(nEntries "$SIMDIR/${J}_sim.root")
  [ "$nsim" -eq "$NEV" ] || { echo "$J GEN_INCOMPLETE ($nsim of $NEV)"; exit 1; }
else
  echo "[$(date +%H:%M:%S)] $J generation already complete ($nsim entries)"
fi

# ---- beam / kinematics truth check, once per sample ---------------------------------------------
# Cheap, and it is the only thing that would catch a swapped target mass or a level that never made
# it into the generator. The (d,d') channel needs the deuteron PDG: with the proton default,
# check 3 finds only stray secondaries and reports a spurious failure.
if ! grep -q "kinematics" "$SIMDIR/${J}_beam.log" 2>/dev/null; then
  root -b -q -l "$REPO/macro/Simulation/ATTPC/17C_dp/check_beam_C17.C(\"$SIMDIR/${J}_sim.root\",2.1290066,$MBEAM,$MTGT,$MRES,$MEJ,$EX,$EJPDG)" > "$SIMDIR/${J}_beam.log" 2>&1 || true
  grep -q "CLOSES" "$SIMDIR/${J}_beam.log" || { echo "$J KINEMATICS_DO_NOT_CLOSE -- see $SIMDIR/${J}_beam.log"; exit 1; }
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
# The fit output must READ BACK with the same number of entries as its input: acceptance_C14.C and
# ex_res_C14_hf.C index the sim and fit trees by entry number and refuse to run if the counts
# disagree, so anything short is not merely incomplete, it is unusable.
if [ "$ngf" -eq "$NEV" ] && grep -q "dE/dx from CATIMA" "$OUT/${JC}_fit.log" 2>/dev/null; then
  echo "[$(date +%H:%M:%S)] $JC fit already complete ($ngf entries)"
else
  echo "[$(date +%H:%M:%S)] $JC genfit at B = $BNEG T, gas $GEO, particle $PART"
  # backwardSeedFix (argument 14) is kFALSE here, matching the 14C(p,p') reference. Unlike the
  # (d,p) arm, the light recoil of an inelastic-scattering channel is FORWARD in the lab
  # (theta_lab = (180 - theta_cm)/2, so 1-85 deg over the whole generated range) and never needs
  # the backward seed reflection.
  root -b -q -l "$UKF/pipeline/fitGenfit_C14.C(\"$JC\",-1,\"$OUT/\",\"\",\"$OUT/\",$BNEG,2,5,\"\",$MEASSIGMA,10.0,170.0,kTRUE,kFALSE,\"$PART\",\"$GEO\",kFALSE,kTRUE,0.0,1,kTRUE,kFALSE,kFALSE,kFALSE)" > "$OUT/${JC}_fit.log" 2>&1
  grep -q "dE/dx from CATIMA" "$OUT/${JC}_fit.log" || { echo "$JC CATIMA_NOT_ENABLED"; exit 1; }
  ngf=$(nEntries "$OUT/${JC}_genfit.root")
  [ "$ngf" -eq "$NEV" ] || { echo "$JC FIT_INCOMPLETE ($ngf of $NEV) -- the file did not close"; exit 1; }
fi

# EJPDG IS NOT OPTIONAL BELOW. acceptance_C14.C and ex_res_C14_hf.C truth-match the light ejectile
# by PDG code and both defaulted to 2212. On (d,d') that matches nothing: genfit fits the deuterons
# perfectly and both macros then report "generated reactions 0 ... acceptance 0.000" on a sample that
# is entirely fine. Both gained an ejPdg argument on 2026-09-02, defaulting to proton so every
# existing (p,p') and (d,p) caller stays byte-identical.
# ---- acceptance -----------------------------------------------------------------------------------
if ! grep -q "overall acceptance" "$OUT/${JC}_acc.log" 2>/dev/null; then
  echo "[$(date +%H:%M:%S)] $JC acceptance"
  root -b -q -l "$PP/acceptance_C14.C(\"$SIMDIR/${J}_sim.root\",\"$OUT/${JC}_genfit.root\",\"$JC\",$EX,$EBEAM,5.0,36,180.0,10.0,0.5,2.0,kTRUE,-1e9,1e9,$MTGT,$MEJ,$MRES,$MBEAM,$EJPDG)" > "$OUT/${JC}_acc.log" 2>&1
  grep -q "overall acceptance" "$OUT/${JC}_acc.log" || { echo "$JC ACC_FAILED"; exit 1; }
  mv -f "$PP/diagnostics/acceptance_${JC}.root" "$OUT/" 2>/dev/null || true
  rm -f "$PP/diagnostics/acceptance_${JC}.png"
fi

# ---- excitation-energy resolution -------------------------------------------------------------------
# ex_res_C14_hf.C hard-codes theta_lab slices 20-90 deg. For THIS channel that is correct and not a
# trap: the light recoil of an inelastic channel cannot exceed 90 deg in the lab, so the window is
# the physical range. (It was wrong for the (d,p) arm, whose transfer peak is backward.) Its flat
# `res` tree carries every accepted event, so inel_summary_C17.C can re-bin without re-running.
if ! grep -q "ex res done" "$OUT/${JC}_exres.log" 2>/dev/null; then
  echo "[$(date +%H:%M:%S)] $JC Ex resolution"
  root -b -q -l "$HF/ex_res_C14_hf.C(\"$SIMDIR/${J}_sim.root\",\"$OUT/${JC}_genfit.root\",\"$JC\",$EX,$EBEAM,5.0,kTRUE,\"$OUT/\",10.0,0.5,2.0,$MTGT,$MEJ,$MRES,$MBEAM,$EJPDG)" > "$OUT/${JC}_exres.log" 2>&1
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
