#!/usr/bin/env bash
# A/B: is CATIMA's per-step dE/dx better than the ASCII table, on 16C(d,t)15C?
#
# WHY THIS CHANNEL. genfit reads the table ONLY below beta*gamma = 0.05, which is KE 3.5 MeV for a
# triton. The (d,t) low branch is 0.8-6 MeV, so these tracks spend most of their path inside the
# region the table serves. The same test on (p,d) is near-blind: those deuterons are ~20 MeV
# against a 2.3 MeV threshold and only dip into it near the endpoint -- measured, 5 of 23 tracks
# moved at all.
#
# THE ONLY DIFFERENCE BETWEEN THE ARMS is the second-to-last argument (catimaELoss). Same reco
# input, same PID gate, same geometry, same par, same measSigma, same backExtrap, same CATIMA MSC
# and straggling in BOTH arms -- so any difference is attributable to the dE/dx source alone.
#
#   arm "tab"  catimaELoss = kFALSE  -> dEdxParam: frozen ASCII curve x ONE GLOBAL density
#   arm "cat"  catimaELoss = kTRUE   -> catima::dedx + dedx_n, using the STEP's own material
#
# A SPREAD OF RUNS, NOT ONE. The first CATIMA recommendation was validated on run_0031 alone and
# looked excellent; run_0031 turned out to be 1 of only 6 runs where matFX worked at all, and the
# other 41 collapsed 40-94%. These 8 span the range, including runs that were formerly the worst
# (0016, 0019, 0098) and formerly the best (0031, 0046).
#
#   ./ab_dt_catimaeloss.sh [nparallel] [nevents]
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
# config.sh reads unset vars: source it BEFORE set -u or it kills the script silently
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export GENFIT=/home/yassid/fair_install/GenFit
export LD_LIBRARY_PATH=$GENFIT/lib:${LD_LIBRARY_PATH:-}
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
set -uo pipefail

NPAR="${1:-4}"
export NEV="${2:--1}"
export ARMS="${ARMS:-tab cat full}"
export RUNS="0016 0019 0031 0040 0046 0057 0084 0098"

export REC=/mnt/f/a1975/reco_d2_dv1104/
export OUT=/mnt/f/a1975/gf_dt_abeloss/
export LOG=/mnt/f/a1975/logs_dt_abeloss/
export PAR=ATTPC.a1975_deuterium_dv1104.par
export GATE=pid/triton_d2_dv1104.json
export TAB=triton_D2_300torr.txt
mkdir -p "$OUT" "$LOG"

# $1 = run number, $2 = arm (tab|cat)
one() {
  n="$1"; arm="$2"; r="run_${n}"
  rc="${REC}${r}_multifit_reco.root"
  [ -s "$rc" ] || { echo "[noreco] $r"; return 0; }
  case "$arm" in
  tab) CEL=kFALSE; CFULL=kFALSE ;;
  cat) CEL=kTRUE;  CFULL=kFALSE ;;
  full) CEL=kTRUE; CFULL=kTRUE ;;
  *) echo "[BUG] unknown arm $arm"; return 1 ;;
  esac
  sfx="_${arm}"
  fo="${OUT}${r}_multifit_genfitter_t${sfx}.root"
  [ -s "${fo}.done" ] && { echo "[have] $r $arm"; return 0; }
  root -l -b -q "fitGenfitter_a1975_deuterium.C(\"${r}_multifit\",${NEV},\"$REC\",\"${sfx}\",\"$OUT\",\
-2.85,2,5,\"$GATE\",4.0,10.0,170.0,kTRUE,kTRUE,1000010030,3.01550072,1,\"t\",\"_reco\",\
\"ATTPC_D300torr_v2_geomanager.root\",kTRUE,0,2,\"$PAR\",kFALSE,kFALSE,kFALSE,\"$TAB\",kTRUE,kTRUE,${CEL},${CFULL})" \
    > "${LOG}gf_${r}${sfx}.log" 2>&1
  if grep -qi 'segmentation violation' "${LOG}gf_${r}${sfx}.log" || [ ! -s "$fo" ]; then
    echo "[FAIL] $r $arm  (see ${LOG}gf_${r}${sfx}.log)"; rm -f "$fo"; return 0
  fi
  # The arms must differ in the dE/dx banner and in NOTHING else. Catch a mislabelled arm here,
  # where it is one line of output, rather than in a spectrum three steps downstream.
  if [ "$arm" = cat ] && ! grep -q "dE/dx from CATIMA" "${LOG}gf_${r}${sfx}.log"; then
    echo "[WARN] $r cat: CATIMA dE/dx banner missing -- arm may be inert"
  fi
  if [ "$arm" = tab ] && grep -q "dE/dx from CATIMA" "${LOG}gf_${r}${sfx}.log"; then
    echo "[WARN] $r tab: CATIMA dE/dx banner PRESENT -- the control arm is contaminated"
  fi
  touch "${fo}.done"; echo "[ok] $r $arm  $(date '+%H:%M:%S')"
}
export -f one

echo "=== (d,t) dE/dx A/B: $(echo $RUNS | wc -w) runs x 2 arms, $NPAR parallel, nev=$NEV ==="
echo "=== control 'tab' = ASCII table; test 'cat' = CATIMA per-step. All else identical. ==="
for a in $ARMS; do for n in $RUNS; do echo "$n $a"; done; done \
  | xargs -P "$NPAR" -I{} bash -c 'one $@' _ {}
echo "=== done: $(ls "$OUT"/*_genfitter_*.root 2>/dev/null | wc -l) fit files in $OUT ==="
