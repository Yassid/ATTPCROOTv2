#!/bin/bash
# Per run: gate (IC 500-900 events + proton-gated tracks) -> UKF fit -> GENFIT fit.
# 4 cores. Outputs local to FITDIR, then rsync to F. Input slim+FRIB local; metadata from F reco.
#   ./fitpipe_C14.sh "run_0055 run_0056 ..." 4
#
# Gas density and geometry are overridable so the same data can be refitted with
# different material. The a1954 gas is 300 torr H2 (rho = 3.553e-5 g/cm3); the original
# production used 6.5e-5 / ATTPC_H600torr, i.e. ~1.9x too much material.
#   DENSITY=3.553e-5 GEONAME=ATTPC_H300torr FITDIR=~/a1954_C14_fit_300torr/ ./fitpipe_C14.sh "..." 4
RUNS="${1:-run_0055}"; NPAR="${2:-4}"
DENSITY="${DENSITY:-3.553e-5}"        # g/cm3, H2 at 300 torr
GEONAME="${GEONAME:-ATTPC_H300torr}"  # genfit navigation/material geometry
REPO="/home/yassid/fair_install/ATTPCROOTv2-OpenKF"
HERE="$REPO/macro/Unpack_HDF5/a1954/UKF"
SLIM="/home/yassid/a1954_C14_reco_hdb_slim/"
FREF="/home/yassid/a1954_C14_reco_hdb/"
FITDIR="${FITDIR:-/home/yassid/a1954_C14_fit/}"; IN="${FITDIR}in/"; LOG="${FITDIR}logs"; mkdir -p "$IN" "$LOG"
echo "[$(date +%H:%M:%S)] fitpipe: density=$DENSITY geom=$GEONAME out=$FITDIR"

source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"

one(){ local r="$1"; local L="$LOG/${r}.log"
  [ -f "$SLIM/${r}_FRIB.root" ] || { echo "no FRIB $r"; return; }
  [ $(stat -c%s "$SLIM/${r}_FRIB.root" 2>/dev/null) -lt 10000 ] && { echo "empty FRIB $r (skip)"; return; }
  # 1) gate
  root -b -q -l "$HERE/pipeline/gate_events_C14.C(\"$r\",\"$SLIM\",\"$IN\",\"$FREF${r}_reco.root\")" > "$L" 2>&1
  [ -f "$IN/${r}_reco.root" ] || { echo "gate failed $r"; return; }
  # 2) UKF
  root -b -q -l "$HERE/pipeline/fitUKF_C14.C(\"$r\",-1,\"proton\",-1,2.85,$DENSITY,\"\",\"$IN\",0.5,0.1,1,10,\"$FITDIR\")" >> "$L" 2>&1
  # 3) GENFIT  (matEffects stays kFALSE here -- change ONE thing at a time; note that with
  #    material effects off the geometry only drives navigation, not energy loss)
  root -b -q -l "$HERE/pipeline/fitGenfit_C14.C(\"$r\",-1,\"$IN\",\"\",\"$FITDIR\",-2.85,2,5,\"\",4.0,10.0,170.0,kFALSE,kFALSE,\"proton\",\"$GEONAME\")" >> "$L" 2>&1
  echo "[$(date +%H:%M:%S)] $r  gated=$(grep -o 'gated-proton events[^,]*' $L|head -1)  ukf=$([ -f $FITDIR${r}_ukf.root ]&&echo ok)  genfit=$([ -f $FITDIR${r}_genfit.root ]&&echo ok)"
}
export -f one; export HERE SLIM FREF FITDIR IN LOG DENSITY GEONAME
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "[$(date +%H:%M:%S)] FIT PIPE DONE. Syncing to F..."
# mirror the local output dir name so a refit never overwrites a previous one on F:
FSYNC="${FSYNC:-/mnt/f/$(basename "${FITDIR%/}")}"
mkdir -p "$FSYNC"
[ -d /mnt/f ] && rsync -a --include='*_ukf.root' --include='*_genfit.root' --exclude='*' "$FITDIR" "$FSYNC/" 2>/dev/null
echo "[$(date +%H:%M:%S)] SYNC DONE -> $FSYNC/"
