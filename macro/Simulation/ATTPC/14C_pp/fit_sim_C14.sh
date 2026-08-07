#!/bin/bash
# UKF- and GENFIT-fit the 14C(p,p') SIMULATION and print the closure result.
#
#   ./fit_sim_C14.sh [NAME] [INDIR] [BSIGN] [BGENFIT]
#     NAME     fit basename, e.g. negBg (default)      INDIR   dir holding <NAME>_reco.root
#     BSIGN    UKF bFieldSign, optional               BGENFIT signed genfit field, optional
#
# FIELD SIGN -- AUTO-DETECTED FROM THE INPUT PATH, because BOTH conventions now exist side by
# side in this tree and a wrong sign does not crash, it silently biases KE:
#     diagnostics/negB/                    B = -28.5 kG (data convention)  ->  -1 / -2.85
#     data/, diagnostics/scan/, ...full/   LEGACY +28.5 kG                 ->  +1 / +2.85
# An unrecognised path is REFUSED unless both signs are given explicitly. The values used and the
# reason are echoed on every run.
#
# NOTHING PRODUCED BEFORE IS INVALIDATED. Every existing sample was generated at +28.5 kG and
# fitted with +1/+2.85 -- internally consistent, and it still reproduces via the LEGACY branch.
# Only what C14_pp_sim.C generates from 2026-08-06 on is in the data convention.
#
# Why the convention was changed: AtSpyralPID::Direction() infers forward/backward purely from
# the sense of rotation and never looks at z, while `polar` comes from d(rho)/dz. At +28.5 kG the
# two disagreed and its consistency check rejected 82 % of sim tracks (10.8 % valid vs 100 % on
# data). Generating in the data convention removes the special case instead of working around it.
set -eo pipefail

NAME=${1:-negBg}
REPO="/home/yassid/fair_install/ATTPCROOTv2-OpenKF"
SIMDIR="$REPO/macro/Simulation/ATTPC/14C_pp"
IN=${2:-"$SIMDIR/diagnostics/negB/"}

# Pick the convention FROM THE SAMPLE, never from a fixed default: the two live side by side
# in this tree and a wrong sign does not crash, it silently biases KE.
#   negB/  -> generated at -28.5 kG (data convention)      -> -1 / -2.85
#   data/, diagnostics/scan/, diagnostics/full/ -> LEGACY +28.5 kG -> +1 / +2.85
case "$IN" in
  *negB*)                              AUTOS="-1"; AUTOG="-2.85"; WHY="B = -28.5 kG, data convention" ;;
  *14C_pp/data/*|*scan/*|*full/*)      AUTOS="+1"; AUTOG="+2.85"; WHY="LEGACY sample generated at +28.5 kG" ;;
  *) echo "REFUSING: cannot tell which field convention $IN was generated with."
     echo "  Pass the sign explicitly: $0 NAME INDIR {-1|+1} {-2.85|+2.85}"
     [ $# -ge 4 ] || exit 2
     AUTOS="$3"; AUTOG="$4"; WHY="caller-specified" ;;
esac
BSIGN=${3:-$AUTOS}
BGENFIT=${4:-$AUTOG}

# --- guard: this script sets a field convention, so refuse anything outside the sim tree.
# It is not a data-fitting script; fitpipe_C14.sh is.
case "$IN" in
  *Simulation/ATTPC/14C_pp/*) : ;;
  *) echo "REFUSING: $IN is outside the 14C_pp simulation tree -- use fitpipe_C14.sh for data"; exit 2 ;;
esac
[ -s "${IN%/}/${NAME}_reco.root" ] || { echo "MISSING ${IN%/}/${NAME}_reco.root"; exit 1; }

set +u # thisroot.sh / config.sh reference unset vars and would trip `set -u`
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
cd "$SIMDIR" || exit 1

UKF="$REPO/macro/Unpack_HDF5/a1954/UKF"
echo "=== fitting $NAME in $IN"
echo "=== UKF bFieldSign=$BSIGN   GENFIT bField=$BGENFIT   [$WHY]"

echo "[$(date +%H:%M:%S)] UKF fit, density 3.553e-5..."
root -b -q -l "$UKF/pipeline/fitUKF_C14.C(\"$NAME\",-1,\"proton\",$BSIGN,2.85,3.553e-5,\"\",\"$IN\",0.5,0.1,1,10,\"$IN\")" 2>&1 | tail -4

echo "[$(date +%H:%M:%S)] GENFIT fit, ATTPC_H300torr, matEffects=kFALSE..."
root -b -q -l "$UKF/pipeline/fitGenfit_C14.C(\"$NAME\",-1,\"$IN\",\"\",\"$IN\",$BGENFIT,2,5,\"\",4.0,10.0,170.0,kFALSE,kFALSE,\"proton\",\"ATTPC_H300torr\")" 2>&1 | tail -3

for F in ukf genfit; do
  echo "[$(date +%H:%M:%S)] excitation spectrum ($F)..."
  root -b -q -l "$UKF/pp/ex_C14.C(\"$NAME\",\"$IN\",161.0,5.0,\"_${NAME}_$F\",1.007825,14.003242,\"\",\"$F\")" 2>&1 | grep "good track"
done

echo "[$(date +%H:%M:%S)] ===== CLOSURE RESULT (truth is exactly 161 MeV -> expect mu~0, slope~0) ====="
for F in ukf genfit; do
  printf "  %-7s " "$F"
  root -b -q -l "$UKF/pp/hyp_C14.C(\"plots/proton_kin_${NAME}_$F.root\",161,\"${NAME}_$F\")" 2>&1 | grep "(0) as measured" | sed 's/(0) as measured *: //'
done
echo "  --- what the sim would ask for to flatten (pure chain bias, truth is known) ---"
for F in ukf genfit; do
  echo "  [$F]"
  root -b -q -l "$UKF/pp/hyp_C14.C(\"plots/proton_kin_${NAME}_$F.root\",161,\"${NAME}_$F\")" 2>&1 | grep -E "^\(A\)|^\(B\)" | sed 's/^/    /'
done
echo "[$(date +%H:%M:%S)] SIM FIT DONE"
