#!/bin/bash
# ============================================================================
# Overnight autonomous 15C+p production (a2091). Launched 2026-07-27, unattended.
#   Phase A  rebuild build_genfit AtUnpack (picks up AtMergerHDFUnpacker)
#   Phase B  reco the 11 raw-FRIBDAQ runs (AtMergerHDFUnpacker, auto-detected)
#   Phase C  ungated UKF + GENFIT(noMat) + GENFIT(matFX) on ALL reco'd runs
# Disk-guarded (stops launching heavy work below MIN_FREE_GB). Resumable
# (skips existing outputs). Fits are UNGATED on purpose — IC/PID gates still
# need a2091 physics tuning; this produces the full three-fitter dataset.
# ============================================================================
set -u
REPO=/home/yassid/fair_install/ATTPCROOTv2
HERE=$REPO/macro/Unpack_HDF5/a2091/UKF
RECO=/home/yassid/a2091_C15_reco
FIT=/home/yassid/a2091_C15_fit
LOG=$RECO/logs
mkdir -p "$FIT" "$LOG"
MASTER=$RECO/overnight.log
MIN_FREE_GB=20
GEOM=ATTPC_H300torr_RT
DENS=3.308e-5

say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$MASTER"; }
free_gb(){ df -BG --output=avail /home | tail -1 | tr -dc '0-9'; }
guard(){ local f; f=$(free_gb); if [ "${f:-0}" -lt "$MIN_FREE_GB" ]; then
           say "DISK GUARD: ${f}GB free < ${MIN_FREE_GB}GB — halting new heavy work."; return 1; fi; return 0; }
srcbuild(){ source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
            source "$REPO/$1/config.sh" >/dev/null 2>&1
            export ROOT_INCLUDE_PATH="$REPO/$1/include:/home/yassid/fair_install/FairRootInstall/include"; }

say "=================== OVERNIGHT 15C+p START (free=$(free_gb)GB) ==================="

# ---- Phase A: build_genfit AtUnpack -------------------------------------------------
say "Phase A: rebuild build_genfit AtUnpack"
( source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
  export FAIRROOTPATH=/home/yassid/fair_install/FairRootInstall SIMPATH=/home/yassid/fair_install/FairSoftInstall
  cmake --build "$REPO/build_genfit" --target AtUnpack -j10 ) > "$LOG/ov_build_genfit_unpack.log" 2>&1
say "Phase A done (exit $?)"
if nm -DC "$REPO/build_genfit/lib/libAtUnpack.so" 2>/dev/null | grep -q AtMergerHDFUnpacker; then
  say "Phase A OK: AtMergerHDFUnpacker present in build_genfit"
else
  say "Phase A WARN: AtMergerHDFUnpacker NOT in build_genfit (fits read reco ROOT, so continuing)"
fi

# ---- Phase B: reco the 11 raw-FRIBDAQ runs -----------------------------------------
# An empty/missing list would make reco_hdb_C15.sh fall back to its run_0138 default,
# which already exists -> the whole phase would no-op silently. Fail loudly instead.
FRIB=$(cat "$HERE/runs_pp_fribdaq.txt" 2>/dev/null)
[ -z "$FRIB" ] && { say "FATAL: runs_pp_fribdaq.txt missing or empty — aborting."; exit 1; }
say "Phase B list: $(echo $FRIB | wc -w) runs"

# Clear the 0-byte reco stubs left by the earlier failed batch so they get redone.
for r in $FRIB; do
  f="$RECO/${r}_reco.root"
  [ -f "$f" ] && [ "$(stat -c%s "$f")" -lt 1000000 ] && { say "rm empty stub $f"; rm -f "$f"; }
done
if guard; then
  say "Phase B: reco FRIBDAQ runs -> $RECO"
  "$HERE/reco_hdb_C15.sh" "$FRIB" 2 >> "$MASTER" 2>&1
  say "Phase B done (free=$(free_gb)GB)"
else
  say "Phase B SKIPPED (disk)"
fi

# ---- Phase C: ungated 3-fitter pass on every reco'd run ----------------------------
ALLRUNS=$(ls "$RECO"/*_reco.root 2>/dev/null | sed 's#.*/##;s/_reco.root//' | sort -u)
say "Phase C: fits over $(echo $ALLRUNS | wc -w) reco'd runs"

# C1: UKF (build/)
srcbuild build
for r in $ALLRUNS; do
  guard || break
  [ "$(stat -c%s "$RECO/${r}_reco.root")" -lt 5000000 ] && continue   # skip junk/corrupt reco
  [ -f "$FIT/${r}_ukf.root" ] && continue
  say "UKF $r"
  root -b -q -l "$HERE/pipeline/fitUKF_C15.C(\"$r\",-1,\"proton\",-1,2.85,$DENS,\"\",\"$RECO/\",0.5,0.1,1,10,\"$FIT/\")" > "$LOG/${r}_ukf.log" 2>&1
done

# C2/C3: GENFIT noMat then matFX (build_genfit)
srcbuild build_genfit
for r in $ALLRUNS; do
  guard || break
  [ "$(stat -c%s "$RECO/${r}_reco.root")" -lt 5000000 ] && continue
  if [ ! -f "$FIT/${r}_genfit_nomat.root" ]; then
    say "GENFIT noMat $r"
    root -b -q -l "$HERE/pipeline/fitGenfit_C15.C(\"$r\",-1,\"$RECO/\",\"_nomat\",\"$FIT/\",-2.85,2,5,\"\",4.0,10.0,170.0,kFALSE,kFALSE,\"proton\",\"$GEOM\")" > "$LOG/${r}_genfit_nomat.log" 2>&1
  fi
  guard || break
  if [ ! -f "$FIT/${r}_genfit_mat.root" ]; then
    say "GENFIT matFX $r"
    root -b -q -l "$HERE/pipeline/fitGenfit_C15.C(\"$r\",-1,\"$RECO/\",\"_mat\",\"$FIT/\",-2.85,2,5,\"\",4.0,10.0,170.0,kTRUE,kFALSE,\"proton\",\"$GEOM\")" > "$LOG/${r}_genfit_mat.log" 2>&1
  fi
done

say "=================== OVERNIGHT 15C+p DONE (free=$(free_gb)GB) ==================="
say "fits in $FIT : ukf=$(ls $FIT/*_ukf.root 2>/dev/null | wc -l)  nomat=$(ls $FIT/*_genfit_nomat.root 2>/dev/null | wc -l)  matFX=$(ls $FIT/*_genfit_mat.root 2>/dev/null | wc -l)"
