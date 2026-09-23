#!/usr/bin/env bash
# genfit the a1975 16C(d,p)17C SIMULATION with the SAME fitter configuration as the adopted DATA
# production (gf_dp_find805), so the resolution it yields describes that production.
#
# Every argument below is copied from prod_dp_findscan.sh's find805 invocation. THREE differ, and
# only three:
#
#   1. *** bField +2.85, NOT -2.85 ***
#      The simulation reverses drift-z against the experiment, so the sense of rotation is
#      opposite. Fitting sim with the data's sign makes anything that infers direction from
#      curvature reject most tracks, and the result then measures the sign error rather than the
#      detector. Consequence to expect downstream: the reconstructed polar angle comes back as
#      180 - true. Do not "fix" that -- it is a property of digitisation.
#      See reference_sim_z_handedness.
#
#   2. *** NO PID GATE ***
#      The data gate (pid/proton_d2_dv1104.json) is drawn on the DATA's Brho/dE-dx plane; the
#      simulation has no ion chamber and its plane is not the same object. The sim is pure g.s.
#      protons by construction, so the gate's job -- species selection -- is already done.
#      BUT NOTE WHAT THIS LETS THROUGH: the reco also finds the BEAM and the 17C RESIDUAL as
#      pattern tracks, and with no gate those get fitted as protons and produce nonsense. They are
#      removed downstream by the kinematic cuts (exsel_dp.C) and can be removed exactly by MC truth
#      (the digitisation ran SetSaveMCInfo, so hits carry trackID/A/Z).
#      fitGenfitter treats an empty gate as "fit everything" SILENTLY -- for data that is a trap
#      (prod_dp_*.sh refuse to start without a gate); here it is the intent, stated out loud.
#
#   3. Input/output paths.
#
# Unchanged from find805: matEffects ON, CATIMA eloss + MSC + straggling ON with
# proton_D2_300torr.txt, measSigma 4.0, theta 10-170, backwardSeedFix ON, pdg 2212,
# ATTPC_D300torr_v2_geomanager.root, and FIND mode 1 with C2MAX 0.1, step 5 %, 80 % -> 5 %,
# minimum 8 clusters.
#
#   ./run_genfit_dp_sim.sh [recoFile]
set -eo pipefail
RECO=${1:-/mnt/f/a1975_C16_dp_sim/gs_s4001_reco_g75000.root}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../../../.." && pwd)
D2="$REPO/macro/Unpack_HDF5/a1975/D2_UKF"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"

[ -s "$RECO" ] || { echo "ERROR: no reco at $RECO"; exit 1; }

IODIR="$(dirname "$RECO")/"
BASE="$(basename "$RECO" .root)"           # gs_s4001_reco_g75000
TAG="${BASE%%_reco_*}"                     # gs_s4001
# NOT `_${BASE#*_}` -- that strips to the FIRST underscore and yields _s4001_reco_g75000, so
# fitGenfitter looked for gs_s4001_s4001_reco_g75000.root. Strip exactly the tag instead.
SUF="${BASE#"$TAG"}"                       # _reco_g75000
OUT="${IODIR}${TAG}_genfitter_p.root"
LOG="${IODIR}${TAG}_genfit.log"

# find805 settings
BFIELD=2.85          # SIM SIGN -- see the header
PAR=ATTPC.a1975_deuterium_dv1104.par
TAB=proton_D2_300torr.txt
GEO=ATTPC_D300torr_v2_geomanager.root
C2MAX=0.1; STEP=5.0; MINPCT=5.0; MAXPCT=80.0; MINCL=8; FINDMODE=1

[ -s "$D2/$TAB" ] || { echo "ERROR: no energy-loss table at $D2/$TAB"; exit 1; }

if [ -f "${OUT}.done" ]; then echo "[$(date +%H:%M:%S)] already done: $OUT"; exit 0; fi

echo "[$(date +%H:%M:%S)] genfit  $BASE   bField=+$BFIELD (SIM)  NO PID GATE  FIND ${MAXPCT}->${MINPCT}%"
echo "                 -> $OUT"
cd "$D2"
# ROOT's exit code is not a verdict -- judge on the readback. See run_gs_dp.sh.
root -l -b -q "fitGenfitter_a1975_deuterium.C(\"$TAG\",-1,\"$IODIR\",\"\",\"$IODIR\",\
$BFIELD,2,5,\"\",4.0,10.0,170.0,kTRUE,kTRUE,2212,1.00782503207,1,\"p\",\"$SUF\",\
\"$GEO\",kTRUE,0,2,\"$PAR\",kFALSE,kFALSE,kFALSE,\"$TAB\",kTRUE,kTRUE,kTRUE,kFALSE,0,kTRUE,\
$C2MAX,kFALSE,$STEP,$MINPCT,$FINDMODE,$MAXPCT,$MINCL)" > "$LOG" 2>&1 || true

if grep -qi 'segmentation violation' "$LOG"; then
  echo "[$(date +%H:%M:%S)] FAIL: segfault during fitting; see $LOG"; rm -f "$OUT"; exit 1
fi
[ -s "$OUT" ] || { echo "[$(date +%H:%M:%S)] FAIL: no output; see $LOG"; exit 1; }
# the two banners that prove the CATIMA backend actually engaged, as the data production checks
grep -q "CATIMA material model: MSC ON, straggling ON" "$LOG" || echo "[WARN] CATIMA MSC banner missing"
grep -q "dE/dx from CATIMA" "$LOG" || echo "[WARN] CATIMA dE/dx banner missing -- fell back to the table"
# `|| true`: the readback root segfaults at exit like every other one here, and under
# `set -eo pipefail` that turns a PASSED verification into a script failure.
vr=$(root -b -l -q "pid/root_ok.C(\"$OUT\")" 2>/dev/null | grep -E '^(VALID|INVALID)' | head -1 || true)
case "$vr" in
  VALID*) touch "${OUT}.done"; echo "[$(date +%H:%M:%S)] DONE  ($vr)  $(du -h "$OUT" | cut -f1)" ;;
  *)      echo "[$(date +%H:%M:%S)] FAIL readback (${vr:-unreadable}) -- removing"; rm -f "$OUT"; exit 1 ;;
esac
