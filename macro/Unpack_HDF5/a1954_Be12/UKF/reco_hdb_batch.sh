#!/bin/bash
# HDBSCAN reco (mover join) for the PID clustering comparison. persistRaw=false.
RUNS="${1:-run_0142}"; NPAR="${2:-2}"
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT="/home/yassid/a1954_Be12_reco_hdb"; LOG="$OUT/logs"; mkdir -p "$LOG"
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh
source /home/yassid/fair_install/ATTPCROOTv2/build/config.sh >/dev/null 2>&1
run_one(){ local r="$1"
  [ -f "$OUT/${r}_reco.root" ] && { echo "skip $r"; return; }
  echo "[$(date +%H:%M:%S)] hdb $r"
  # unpackReco_Be12(file,nEv,persistRaw,outDir,filepath,doSC,applyTimeCorr,psaType,primSigma,thr,praType)
  root -b -q -l "$HERE/pipeline/unpackReco_Be12.C(\"$r\",0,false,\"$OUT/\",\"/media/yassid/NSCL_e15250/data/a1954_remerged/h5/\",false,false,\"multifit\",0,20,\"hdbscan\")" > "$LOG/${r}_reco.log" 2>&1
  echo "[$(date +%H:%M:%S)] done $r exit $? ($(ls -la $OUT/${r}_reco.root 2>/dev/null|awk '{printf "%d MB",$5/1e6}'))"
}
export -f run_one; export OUT HERE LOG
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'run_one "$@"' _ {}
echo "HDBSCAN reco batch done."
