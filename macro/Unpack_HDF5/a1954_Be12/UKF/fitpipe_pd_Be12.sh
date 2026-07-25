#!/bin/bash
# 12Be(p,d)11Be pass: per run  gate (IC 500-900 + DEUTERON-gated tracks) -> UKF(deuteron).
# Separate dirs from the (p,p') pass so nothing is overwritten. GENFIT is skipped on
# purpose (UKF beats it on this data set, see ANALYSIS_REPORT).
#   ./fitpipe_pd_Be12.sh "run_0143 run_0147" 4 [thMin]
RUNS="${1:-run_0143}"; NPAR="${2:-4}"; THMIN="${3:-0}"
REPO="/home/yassid/fair_install/ATTPCROOTv2-OpenKF"
HERE="$REPO/macro/Unpack_HDF5/a1954_Be12/UKF"
SLIM="/home/yassid/a1954_Be12_reco_hdb_slim/"
FREF="/mnt/f/a1954_Be12_reco_hdb/"
FITDIR="/home/yassid/a1954_Be12_fit_pd/"; IN="${FITDIR}in/"; LOG="${FITDIR}logs"; mkdir -p "$IN" "$LOG"
DGATE="$HERE/pid/deuteron_12Be.json"

source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"

one(){ local r="$1"; local L="$LOG/${r}.log"
  [ -f "$SLIM/${r}_FRIB.root" ] || { echo "no FRIB $r"; return; }
  [ $(stat -c%s "$SLIM/${r}_FRIB.root" 2>/dev/null) -lt 10000 ] && { echo "empty FRIB $r (skip)"; return; }
  # 1) gate: IC window + deuteron PID polygon, theta_lab > THMIN
  root -b -q -l "$HERE/pipeline/gate_events_Be12.C(\"$r\",\"$SLIM\",\"$IN\",\"$FREF${r}_reco.root\",500,900,1050,1250,2.85,\"$DGATE\",$THMIN)" > "$L" 2>&1
  [ -f "$IN/${r}_reco.root" ] || { echo "gate failed $r"; return; }
  # 2) UKF with the DEUTERON mass hypothesis
  root -b -q -l "$HERE/pipeline/fitUKF_Be12.C(\"$r\",-1,\"deuteron\",-1,2.85,6.5e-5,\"\",\"$IN\",0.5,0.1,1,10,\"$FITDIR\")" >> "$L" 2>&1
  echo "[$(date +%H:%M:%S)] $r  $(grep -o '[0-9]* gated-proton events, [0-9]* tracks' $L|head -1)  ukf=$([ -f $FITDIR${r}_ukf.root ]&&echo ok||echo FAIL)"
}
export -f one; export HERE SLIM FREF FITDIR IN LOG DGATE THMIN
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "[$(date +%H:%M:%S)] (p,d) FIT PIPE DONE -> $FITDIR"
