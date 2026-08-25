#!/bin/bash
# 14C(p,p') GENFIT refit with MATERIAL EFFECTS + native CATIMA dE/dx, 2026-08-25.
#
# What changes vs the adopted gf_xtr production (26 Jul, matEffects=kFALSE + a manual CATIMA
# term over the vertex gap):
#   * matEffects = kTRUE, so genfit models the gas along the WHOLE track, not just the gap;
#   * dE/dx from CATIMA per step and per material (catimaELoss), which is the a1975 adopted
#     setting -- it more than halved the collapsed-fit rate there at unchanged energy scale;
#   * geometry ATTPC_H300torr_RT (3.308e-5), NOT ATTPC_H300torr (3.553e-5, a 0 C value):
#     with matEffects on the geometry IS the material model, so the 7.4 % matters;
#   * matFallback = kFALSE, so a track that throws is DROPPED rather than silently refitted
#     without material and kept -- otherwise one spectrum mixes two physics models;
#   * NO manual eloss: with matEffects on, extrapolateToLine already integrates the vertex gap
#     and the manual term would count it twice (the macro now refuses it outright).
#
# ONLY SAFE AFTER 2026-08-17 c3e191ce (TGeo navigator re-seed). Before it, matFX-on collapsed
# 46 % of the 14C fits -- measured on a1954_C14_fit_300torr_matfx.
#
#   ./genfit_catima_C14.sh "run_0055 run_0056 ..." 4
RUNS="${1:-run_0055}"; NPAR="${2:-4}"
REPO="/home/yassid/fair_install/ATTPCROOTv2-OpenKF"
HERE="$REPO/macro/Unpack_HDF5/a1954/UKF"
IN="${IN:-/home/yassid/a1954_C14_fit_300torr/in/}"
OUTDIR="${OUTDIR:-/home/yassid/a1954_C14_gf_catima/}"; LOG="$OUTDIR/logs"; mkdir -p "$LOG"
GEONAME="${GEONAME:-ATTPC_H300torr_RT}"

source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
echo "[$(date +%H:%M:%S)] genfit matFX+CATIMA  geom=$GEONAME  in=$IN  out=$OUTDIR"

one(){ local r="$1"; local L="$LOG/${r}.log"
  [ -f "$IN/${r}_reco.root" ] || { echo "no gated input $r"; return; }
  # args: file,nEv,ioDir,suffix,outDir,B,minIt,maxIt,pidGate,measSigma,thMin,thMax,
  #       matEffects,backwardSeedFix,particle,geoName,matFallback,backExtrap,manualEloss,matA,
  #       catimaELoss,catimaELossFull,catimaMSC,catimaStraggling
  root -b -q -l "$HERE/pipeline/fitGenfit_C14.C(\"$r\",-1,\"$IN\",\"\",\"$OUTDIR\",-2.85,2,5,\"\",4.0,10.0,170.0,kTRUE,kFALSE,\"proton\",\"$GEONAME\",kFALSE,kTRUE,0.0,1,kTRUE,kFALSE,kFALSE,kFALSE)" > "$L" 2>&1
  echo "[$(date +%H:%M:%S)] $r genfit=$([ -f $OUTDIR${r}_genfit.root ]&&echo ok||echo FAIL)"
}
export -f one; export HERE IN OUTDIR LOG GEONAME
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "[$(date +%H:%M:%S)] DONE -> $OUTDIR"
