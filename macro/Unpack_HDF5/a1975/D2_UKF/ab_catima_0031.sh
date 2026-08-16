#!/usr/bin/env bash
# Material-model arms for the a1975 (d,t) genfit fit, run on a single run so they can be
# differenced track-by-track. Each arm differs from the others in ONE setting; everything else
# -- gate, geometry, beam, gas density, iteration counts, seed handling -- is identical.
#
# THREE SETTINGS THAT ARE NOT THE PRODUCTION DEFAULTS, and why:
#
#   matEffects = kTRUE for every arm except `nomat`.  fit_dt_dv1104.sh runs kFALSE.
#       noiseCoulomb and noiseBetheBloch are the ONLY places the CATIMA backend is consulted and
#       they are never called with material effects off, so an A/B on the production default
#       would return a perfect null and read as "the backend does nothing" rather than "the
#       backend was never reached".
#
#   matFallback = kFALSE.  AtGenfitter defaults this to true, silently refitting a track whose
#       material-effects fit threw with setNoEffects(true) and keeping it. A production built
#       that way is a BLEND of two physics models, which inflates the width and can fake a
#       shift -- exactly what must not happen in a comparison whose whole subject is the
#       material model. With it off, every surviving track is pure material-effects, and the
#       yields become a second observable.
#
#   eLossTable held FIXED across the arms that are being differenced.  Below beta*gamma = 0.05
#       (KE = 3.5 MeV for a triton) genfit applies NO stopping power without a curve loaded, and
#       the (d,t) low branch is 0.8-6 MeV. Held constant it cannot confound; the `notable` arm
#       varies it deliberately, and only against `on`.
#
# THE ARMS
#   nomat    matFX off                             -- what production runs today
#   off      matFX on, genfit Highland             -- collapses 62% of fits, see the cmp macro
#   on       matFX on, CATIMA MSC + straggling     -- the candidate
#   nostrag  matFX on, CATIMA MSC only             -- isolates the straggling term
#   notable  matFX on, CATIMA MSC + straggling, NO dE/dx table -- isolates the table
#
# There is deliberately NO "straggling only" arm: turning CATIMA MSC off reverts to genfit's
# Highland, which collapses ndf on 62% of tracks, so such an arm cannot produce a clean sample
# to difference against. Straggling is the separable term, MSC is not.
#
#   ./ab_catima_0031.sh [run] [arm ...]        default: all five
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
# config.sh reads unset vars: source it BEFORE set -u or it kills the script silently
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export GENFIT=/home/yassid/fair_install/GenFit
export LD_LIBRARY_PATH=$GENFIT/lib:${LD_LIBRARY_PATH:-}
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
set -uo pipefail

RUN="${1:-0031}"; shift || true
ARMS="${*:-nomat off on nostrag notable}"
REC=/mnt/f/a1975/reco_d2_dv1104/
PAR=ATTPC.a1975_deuterium_dv1104.par
GATE=pid/triton_d2_dv1104.json
RHO=6.61e-5
TAB=triton_D2_300torr.txt
LOG=/mnt/f/a1975/logs_ab_catima/
mkdir -p "$LOG"

# name -> matFX  catimaMSC  catimaStraggling  table
spec() {
  case "$1" in
    nomat)   echo "kFALSE kFALSE kFALSE $TAB" ;;
    off)     echo "kTRUE  kFALSE kFALSE $TAB" ;;
    on)      echo "kTRUE  kTRUE  kTRUE  $TAB" ;;
    nostrag) echo "kTRUE  kTRUE  kFALSE $TAB" ;;
    notable) echo "kTRUE  kTRUE  kTRUE  ''"   ;;
    *)       echo "" ;;
  esac
}

arm() {
  tag="$1"
  read -r MFX MSC STR TB <<<"$(spec "$tag")"
  [ -n "${MFX:-}" ] || { echo "[skip] unknown arm $tag"; return 0; }
  [ "$TB" = "''" ] && TB=""
  out="/mnt/f/a1975/gf_dt_ab_${tag}/"
  mkdir -p "$out"
  f="${out}run_${RUN}_multifit_genfitter_t.root"
  [ -s "$f" ] && { echo "[have] $tag"; return 0; }
  echo "[$(date '+%H:%M:%S')] start $tag  (matFX=$MFX msc=$MSC strag=$STR table='${TB}')"
  root -l -b -q "fitGenfitter_a1975_deuterium.C(\"run_${RUN}_multifit\",-1,\"$REC\",\"\",\"$out\",\
-2.85,2,5,\"$GATE\",4.0,10.0,170.0,${MFX},kTRUE,1000010030,3.01550072,1,\"t\",\"_reco\",\
\"ATTPC_D300torr_v2_geomanager.root\",kTRUE,${RHO},2,\"$PAR\",kFALSE,kFALSE,kFALSE,\"$TB\",${MSC},${STR})" \
    > "${LOG}${tag}_run_${RUN}.log" 2>&1
  if grep -qi 'segmentation violation' "${LOG}${tag}_run_${RUN}.log" || [ ! -s "$f" ]; then
    echo "[$(date '+%H:%M:%S')] FAIL $tag  (see ${LOG}${tag}_run_${RUN}.log)"
  else
    echo "[$(date '+%H:%M:%S')] done $tag  $(du -h "$f" | cut -f1)"
  fi
}

for a in $ARMS; do arm "$a" & done
wait
echo "=== arms finished: $ARMS ==="
