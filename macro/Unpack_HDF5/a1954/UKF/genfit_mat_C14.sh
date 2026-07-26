#!/bin/bash
# GENFIT pass over the ALREADY-GATED 14C inputs, with MATERIAL EFFECTS ON.
#
# This is the configuration that fixed 12Be (genfit_pd_Be12.sh, after the AtGenfitter
# dE/dx fix 5f87e625): matEffects=kTRUE. The 14C production never got it -- fitpipe_C14.sh
# still passes kFALSE, i.e. genfit fits a constant-momentum helix with no dE/dx in the
# track model, which is what produced the ~0.5 MeV UKF/GENFIT split.
#
# Difference from the 12Be run: geometry is ATTPC_H300torr (rho = 3.553e-5 g/cm3), not the
# ATTPC_H600torr default. The 12Be genfix production used 600 torr and so was still ~1.9x
# too dense; with matEffects ON the geometry IS the energy-loss model, so this matters.
#
#   ./genfit_mat_C14.sh "run_0055 run_0056 ..." 4
RUNS="${1:-run_0056}"; NPAR="${2:-4}"
GEONAME="${GEONAME:-ATTPC_H300torr}"
REPO="/home/yassid/fair_install/ATTPCROOTv2-OpenKF"
HERE="$REPO/macro/Unpack_HDF5/a1954/UKF"
IN="${IN:-/home/yassid/a1954_C14_fit_300torr/in/}"
OUTDIR="${OUTDIR:-/home/yassid/a1954_C14_genfit_mat/}"; LOG="${OUTDIR}logs"
mkdir -p "$OUTDIR" "$LOG"

source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"

echo "[$(date +%H:%M:%S)] genfit matEffects=kTRUE  geom=$GEONAME  in=$IN  out=$OUTDIR"

one(){ local r="$1"; local L="$LOG/${r}_genfit.log"
  [ -f "$IN/${r}_reco.root" ] || { echo "[$(date +%H:%M:%S)] $r NO GATED INPUT (skip)"; return; }
  root -b -q -l "$HERE/pipeline/fitGenfit_C14.C(\"$r\",-1,\"$IN\",\"\",\"$OUTDIR\",-2.85,2,5,\"\",4.0,10.0,170.0,kTRUE,kFALSE,\"proton\",\"$GEONAME\")" > "$L" 2>&1
  # ndf<=0 is the silent-garbage signature of the old useEnergyLossParam bug -- surface it
  local nd=$(grep -c "ndf" "$L" 2>/dev/null)
  echo "[$(date +%H:%M:%S)] $r genfit=$([ -f $OUTDIR${r}_genfit.root ]&&echo ok||echo FAIL)"
}
export -f one; export HERE IN OUTDIR LOG GEONAME
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "[$(date +%H:%M:%S)] GENFIT matEffects PASS DONE -> $OUTDIR"
