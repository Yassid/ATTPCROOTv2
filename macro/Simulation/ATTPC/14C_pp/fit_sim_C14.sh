#!/bin/bash
# Wait for run_reco_C14.C to finish, then UKF-fit the SIMULATION and print the closure result.
#
# bFieldSign = +1 -- SIMULATION ONLY. Experimental data uses -1 (the fitUKF_C14.C default,
# and what fitpipe_C14.sh passes). The sim reverses drift-z handedness in digitization and
# the experiment does not. Passing +1 on real data does not crash, it silently biases KE.
# This script therefore refuses to run on anything that is not the sim file.
set -u
REPO="/home/yassid/fair_install/ATTPCROOTv2-OpenKF"
SIMDIR="$REPO/macro/Simulation/ATTPC/14C_pp"
RECOLOG="/home/yassid/c14sim_reco.log"
IN="$SIMDIR/data/"
NAME="sim"
BSIGN="+1"

# --- guard: refuse to apply the simulation handedness to anything but the sim ----------
case "$IN" in
  *Simulation/ATTPC/14C_pp/data/) : ;;
  *) echo "REFUSING: bFieldSign=$BSIGN is simulation-only, but input dir is $IN"; exit 2 ;;
esac
if [ "$NAME" != "sim" ]; then echo "REFUSING: $NAME does not look like the simulation"; exit 2; fi

echo "[$(date +%H:%M:%S)] waiting for the sim reco..."
while true; do
  grep -q "sim reco done" "$RECOLOG" 2>/dev/null && break
  pgrep -f "run_reco_C14.C" >/dev/null 2>&1 || { echo "reco exited without completing"; exit 1; }
  sleep 20
done
echo "[$(date +%H:%M:%S)] reco done: $(grep -oE 'Real [0-9.]+ s' "$RECOLOG" | tail -1)"

set +u # thisroot.sh / config.sh reference unset vars and would trip `set -u`
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
cd "$SIMDIR" || exit 1

UKF="$REPO/macro/Unpack_HDF5/a1954/UKF"
echo "[$(date +%H:%M:%S)] UKF fit (bFieldSign=$BSIGN, SIMULATION), density 3.553e-5..."
root -b -q -l "$UKF/pipeline/fitUKF_C14.C(\"$NAME\",-1,\"proton\",$BSIGN,2.85,3.553e-5,\"\",\"$IN\",0.5,0.1,1,10,\"$IN\")" 2>&1 | tail -4

# GENFIT on the same tracks. bField is signed directly here (not a +-1 flag): the data uses
# -2.85, so the simulation handedness is +2.85. matEffects=kFALSE and ATTPC_H300torr, i.e.
# exactly the production configuration.
echo "[$(date +%H:%M:%S)] GENFIT fit (bField=+2.85, SIMULATION), ATTPC_H300torr, matEffects=kFALSE..."
root -b -q -l "$UKF/pipeline/fitGenfit_C14.C(\"$NAME\",-1,\"$IN\",\"\",\"$IN\",+2.85,2,5,\"\",4.0,10.0,170.0,kFALSE,kFALSE,\"proton\",\"ATTPC_H300torr\")" 2>&1 | tail -3

for F in ukf genfit; do
  echo "[$(date +%H:%M:%S)] excitation spectrum ($F)..."
  root -b -q -l "$UKF/pp/ex_C14.C(\"$NAME\",\"$IN\",161.0,5.0,\"_sim_$F\",1.007825,14.003242,\"\",\"$F\")" 2>&1 | grep "good track"
done

echo "[$(date +%H:%M:%S)] ===== CLOSURE RESULT (truth is exactly 161 MeV -> expect mu~0, slope~0) ====="
for F in ukf genfit; do
  printf "  %-7s " "$F"
  root -b -q -l "$UKF/pp/hyp_C14.C(\"plots/proton_kin_sim_$F.root\",161,\"sim_$F\")" 2>&1 | grep "(0) as measured" | sed 's/(0) as measured *: //'
done
echo "  --- what the sim would ask for to flatten (pure chain bias, truth is known) ---"
for F in ukf genfit; do
  echo "  [$F]"
  root -b -q -l "$UKF/pp/hyp_C14.C(\"plots/proton_kin_sim_$F.root\",161,\"sim_$F\")" 2>&1 | grep -E "^\(A\)|^\(B\)" | sed 's/^/    /'
done
echo "[$(date +%H:%M:%S)] SIM FIT DONE"
