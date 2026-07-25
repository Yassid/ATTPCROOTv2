#!/bin/bash
# HDBSCAN reco for a1954 12Be, reading AND writing on drive F (WD_Cris).
# Config: AtPSAMultiFit (thr20) + AtDirDeDxCleaner + AtTrackFinderHDBSCAN (mover), persistRaw=false.
# Corrected env vs reco_hdb_batch.sh: this repo's build + canonical build/include ROOT_INCLUDE_PATH.
#
#   ./reco_hdb_F.sh "run_0142 run_0143 ..." 2
RUNS="${1:-run_0142}"; NPAR="${2:-2}"
REPO="/home/yassid/fair_install/ATTPCROOTv2-OpenKF"
HERE="$REPO/macro/Unpack_HDF5/a1954_Be12/UKF"
IN="/mnt/f/a1954_remerged/"
OUT="/mnt/f/a1954_Be12_reco_hdb"; LOG="$OUT/logs"; mkdir -p "$LOG"

source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"

run_one(){ local r="$1"
  [ -f "$OUT/${r}_reco.root" ] && { echo "skip $r (exists)"; return; }
  [ -f "$IN/${r}.h5" ] || { echo "MISSING input $r"; return; }
  echo "[$(date +%H:%M:%S)] start $r"
  # unpackReco_Be12(file,nEv,persistRaw,outDir,filepath,doSC,applyTimeCorr,psaType,primSigma,thr,praType)
  root -b -q -l "$HERE/pipeline/unpackReco_Be12.C(\"$r\",0,false,\"$OUT/\",\"$IN\",false,false,\"multifit\",0,20,\"hdbscan\")" > "$LOG/${r}_reco.log" 2>&1
  local ec=$?
  echo "[$(date +%H:%M:%S)] done $r exit $ec ($(ls -la $OUT/${r}_reco.root 2>/dev/null | awk '{printf "%d MB",$5/1e6}'))"
}
export -f run_one; export OUT IN HERE LOG
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'run_one "$@"' _ {}
echo "[$(date +%H:%M:%S)] HDBSCAN reco batch DONE."
