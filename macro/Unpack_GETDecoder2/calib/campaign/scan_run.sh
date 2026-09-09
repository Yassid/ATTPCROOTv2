#!/bin/bash
# Dump one unpacked run to a slim CSV, select alpha+alpha elastic candidates, keep the pickle.
#   scan_run.sh <run> <setting>
RUN=$1; SET=$2
BASE=/home/yassid/dec2014_calib/campaign
OUT="/media/yassid/Seagate Hub/ATTPC/Dec2014_alphas"
ROOTDIR=/home/yassid/fair_install/ATTPCROOTv2_fr19port
SCRATCH=$BASE/tmp; mkdir -p "$SCRATCH" "$BASE/$SET"
PKL="$BASE/$SET/cands_$RUN.pkl"
[ -s "$PKL" ] && { echo "[skip] $RUN"; exit 0; }
SRC="$OUT/$SET/alpha_run_${RUN}_hits.root"
[ -s "$SRC" ] || { echo "[miss] $RUN no unpacked file"; exit 1; }

# one thread per worker: the driver already runs many runs in parallel
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
source $ROOTDIR/setup_fr19port.sh >/dev/null 2>&1
cd $ROOTDIR/macro/Unpack_GETDecoder2 || exit 1
FULL=$SCRATCH/full_$RUN.csv; SLIM=$SCRATCH/slim_$RUN.csv
root -l -b -q "dump_hits.C(\"$SRC\",\"$FULL\",\"hits\",-1)" > $SCRATCH/dump_$RUN.log 2>&1
awk -F, 'NR==1{print "event,tb,x,y,q"; next}{print $1","$3","$4","$5","$10}' "$FULL" > "$SLIM"
rm -f "$FULL"
python3 $BASE/select_run.py "$SLIM" "$PKL" 2>>$SCRATCH/sel_$RUN.log
rm -f "$SLIM"
