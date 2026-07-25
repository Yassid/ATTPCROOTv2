#!/bin/bash
# Per run: gate (IC 500-900 events + proton-gated tracks) -> UKF fit -> GENFIT fit.
# 4 cores. Outputs local to FITDIR, then rsync to F. Input slim+FRIB local; metadata from F reco.
#   ./fitpipe_Be12.sh "run_0142 run_0143 ..." 4
RUNS="${1:-run_0142}"; NPAR="${2:-4}"
REPO="/home/yassid/fair_install/ATTPCROOTv2-OpenKF"
HERE="$REPO/macro/Unpack_HDF5/a1954_Be12/UKF"
SLIM="/home/yassid/a1954_Be12_reco_hdb_slim/"
FREF="/mnt/f/a1954_Be12_reco_hdb/"
FITDIR="/home/yassid/a1954_Be12_fit/"; IN="${FITDIR}in/"; LOG="${FITDIR}logs"; mkdir -p "$IN" "$LOG"

source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"

one(){ local r="$1"; local L="$LOG/${r}.log"
  [ -f "$SLIM/${r}_FRIB.root" ] || { echo "no FRIB $r"; return; }
  [ $(stat -c%s "$SLIM/${r}_FRIB.root" 2>/dev/null) -lt 10000 ] && { echo "empty FRIB $r (skip)"; return; }
  # 1) gate
  root -b -q -l "$HERE/pipeline/gate_events_Be12.C(\"$r\",\"$SLIM\",\"$IN\",\"$FREF${r}_reco.root\")" > "$L" 2>&1
  [ -f "$IN/${r}_reco.root" ] || { echo "gate failed $r"; return; }
  # 2) UKF
  root -b -q -l "$HERE/pipeline/fitUKF_Be12.C(\"$r\",-1,\"proton\",-1,2.85,6.5e-5,\"\",\"$IN\",0.5,0.1,1,10,\"$FITDIR\")" >> "$L" 2>&1
  # 3) GENFIT
  root -b -q -l "$HERE/pipeline/fitGenfit_Be12.C(\"$r\",-1,\"$IN\",\"\",\"$FITDIR\",-2.85,2,5,\"\",4.0,10.0,170.0,kFALSE,kFALSE,\"proton\")" >> "$L" 2>&1
  echo "[$(date +%H:%M:%S)] $r  gated=$(grep -o 'gated-proton events[^,]*' $L|head -1)  ukf=$([ -f $FITDIR${r}_ukf.root ]&&echo ok)  genfit=$([ -f $FITDIR${r}_genfit.root ]&&echo ok)"
}
export -f one; export HERE SLIM FREF FITDIR IN LOG
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "[$(date +%H:%M:%S)] FIT PIPE DONE. Syncing to F..."
mkdir -p /mnt/f/a1954_Be12_fit
rsync -a --include='*_ukf.root' --include='*_genfit.root' --exclude='*' "$FITDIR" /mnt/f/a1954_Be12_fit/ 2>/dev/null
echo "[$(date +%H:%M:%S)] SYNC DONE -> /mnt/f/a1954_Be12_fit/"
