#!/bin/bash
# Unpack the FRIB/IC group from the a1954 HDF5 on F -> local <run>_FRIB.root (small).
#   ./frib_F_batch.sh "run_0142 run_0143 ..." 4
RUNS="${1:-run_0142}"; NPAR="${2:-4}"
REPO="/home/yassid/fair_install/ATTPCROOTv2-OpenKF"
HERE="$REPO/macro/Unpack_HDF5/a1954_Be12/UKF"
IN="/mnt/f/a1954_remerged/"
OUT="/home/yassid/a1954_Be12_reco_hdb_slim/"; LOG="${OUT}logs"; mkdir -p "$LOG"

source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"

frib_one(){ local r="$1"
  [ -f "$OUT/${r}_FRIB.root" ] && { echo "skip $r"; return; }
  [ -f "$IN/${r}.h5" ] || { echo "no input $r"; return; }
  root -b -q -l "$HERE/pipeline/unpackFRIB_Be12.C(\"$r\",\"$OUT\",\"$IN\")" > "$LOG/${r}_FRIB.log" 2>&1
  echo "[$(date +%H:%M:%S)] $r exit $? ($(ls -la $OUT/${r}_FRIB.root 2>/dev/null | awk '{printf "%d MB",$5/1e6}'))"
}
export -f frib_one; export OUT IN HERE LOG
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'frib_one "$@"' _ {}
echo "[$(date +%H:%M:%S)] FRIB unpack DONE."
