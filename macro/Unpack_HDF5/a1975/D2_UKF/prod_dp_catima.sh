#!/usr/bin/env bash
# a1975 16C(d,p)17C genfit production at dv 1.10424, matEffects ON, CATIMA backend.
#
# This is prod_dt_catima.sh with the SPECIES changed and NOTHING ELSE. That is deliberate: the
# (d,t) header justifies every flag on measurements taken on these same 47 D2 runs, and this
# channel shares the reco, the geometry, the par and the ion chamber with it. Only the four
# species arguments, the gate, the dE/dx table and the output paths differ.
#
#   pdg 1000010030 -> 2212, massAmu 3.01550072 -> 1.00782503207, speciesTag "t" -> "p"
#       These are also fitGenfitter_a1975_deuterium.C's own DEFAULTS, which is exactly the trap
#       that killed the previous (d,p) production: overnight_full_dp.sh passed only 6 arguments,
#       so it silently took the defaults for everything after them -- matEffects kFALSE, geoName
#       ATTPC_H1bar (ordinary hydrogen, not D2), the 1.136-baseline par, matFallback kTRUE.
#       The species happened to be right; the physics was not. Every argument is passed here.
#
#   GATE  pid/proton_d2_dv1104.json, drawn on the species-blind plane pid_plane_dt_dv1104.root
#       (the "dt" in that filename is a misnomer -- pid_plane_dt.C runs AtSpyralPID::Estimate
#       with no species hypothesis and no cuts, so it is the plane for every D2 species).
#       LOOSE by policy: genfit gates BEFORE fitting, so whatever runs here is frozen into the
#       output. Fitting the superset keeps the tight gate re-cuttable at the cache stage, where
#       brho / dedx / sqrtdedx / icmax are stored per track.
#
#   TAB   proton_D2_300torr.txt, generated 2026-09-05 by
#         make_eloss_table.C("proton_D2_300torr.txt", 1, 1, 1.00782503207, 6.61e-5, 2.014)
#       Below beta*gamma = 0.05 genfit applies NO stopping power without a curve loaded. For a
#       PROTON that threshold is KE = 1.17 MeV, not the triton's 3.5, so the dead region is much
#       smaller here -- but backward (d,p) protons land in it, and those are the transfer tracks.
#       VALIDATED against the production triton table: at equal velocity (KE_t = 2.994 x KE_p)
#       the two agree to 0.04% over 0.2-10 MeV, which is what Z^2 scaling demands for equal Z.
#
#   OUT   gf_dp_cateloss/ + logs_dp_cateloss/, alongside the (d,t) products. Nothing is
#       overwritten; the deleted 2026-06-30 production lived in reco_d2/ and does not come back.
#
# UNCHANGED and NOT to be "fixed" without a measurement: bField -2.85, iterations 2/5,
# measSigma 4.0, theta 10-170 (protons go BOTH ways, so backwardSeedFix kTRUE matters here more
# than it does for tritons), geometry ATTPC_D300torr_v2, par ..._dv1104, manualElossDensity 0
# (extrapolateToLine already integrates the vertex gap under matEffects -- passing 6.61e-5 double
# counts it), matFallback kFALSE, catimaELoss kTRUE / catimaELossFull kFALSE.
#
#   ./prod_dp_catima.sh [nparallel] [runs...]
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
# config.sh reads unset vars: source it BEFORE set -u or it kills the script silently
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export GENFIT=/home/yassid/fair_install/GenFit
export LD_LIBRARY_PATH=$GENFIT/lib:${LD_LIBRARY_PATH:-}
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
set -uo pipefail

NPAR="${1:-6}"; shift || true
RUNS="${*:-0016 0017 0018 0019 0020 0021 0022 0023 0026 0027 0031 0032 0034 0036 0037 0038 0039 0040 0041 0042 0043 0044 0046 0048 0057 0058 0076 0077 0078 0079 0080 0082 0083 0084 0085 0086 0087 0088 0089 0091 0092 0095 0096 0097 0098 0102 0103}"

export REC=/mnt/f/a1975/reco_d2_dv1104/
export GF="${GF:-/mnt/f/a1975/gf_dp_cateloss/}"
export LOG="${LOG:-/mnt/f/a1975/logs_dp_cateloss/}"
export PAR=ATTPC.a1975_deuterium_dv1104.par
export GATE="${GATE:-pid/proton_d2_dv1104.json}"
export TAB=proton_D2_300torr.txt
export CEL="${CEL:-kTRUE}"
export CFULL="${CFULL:-kFALSE}"

# A missing gate would NOT stop the fitter -- fitGenfitter treats an empty/absent pidGate as
# "no gate" and happily fits every track in the file, which is a 47-run production of garbage
# that looks successful. Fail loudly instead.
[ -s "$GATE" ] || { echo "ERROR: PID gate '$GATE' missing or empty. Draw it with"; \
                    echo "         root -l 'pid/run_gate_dp.C'"; exit 1; }
[ -s "$TAB" ]  || { echo "ERROR: dE/dx table '$TAB' missing."; exit 1; }
mkdir -p "$GF" "$LOG"

one() {
  n="$1"; r="run_${n}"
  rc="${REC}${r}_multifit_reco.root"
  [ -s "$rc" ] || { echo "[noreco] $r"; return 0; }
  fo="${GF}${r}_multifit_genfitter_p.root"
  # RESUME GUARD. It was `[ -s "${fo}.done" ]`, inherited from prod_dt_catima.sh, and it could
  # NEVER fire: the marker is made with `touch`, so it is zero bytes, and -s tests NON-EMPTY.
  # Every rerun therefore refitted all 47 runs from scratch while printing nothing to say so.
  # -e fixes that, but the marker only records that the job returned -- it says nothing about
  # whether the file is usable, and a killed job leaves a header-sized stub that -e is happy with.
  # So read the product back, which is the only test that distinguishes the two.
  if [ -e "${fo}.done" ] && [ -s "$fo" ]; then
    vr=$(root -b -l -q "pid/root_ok.C(\"$fo\")" 2>/dev/null | grep -E '^(VALID|INVALID)' | head -1)
    case "$vr" in
      VALID*) echo "[have] $r  (${vr})"; return 0 ;;
      *)      echo "[redo] $r  (marker present but file is ${vr:-unreadable})"; rm -f "${fo}.done" ;;
    esac
  fi
  root -l -b -q "fitGenfitter_a1975_deuterium.C(\"${r}_multifit\",-1,\"$REC\",\"\",\"$GF\",\
-2.85,2,5,\"$GATE\",4.0,10.0,170.0,kTRUE,kTRUE,2212,1.00782503207,1,\"p\",\"_reco\",\
\"ATTPC_D300torr_v2_geomanager.root\",kTRUE,0,2,\"$PAR\",kFALSE,kFALSE,kFALSE,\"$TAB\",kTRUE,kTRUE,${CEL},${CFULL})" \
    > "${LOG}gf_${r}.log" 2>&1
  if grep -qi 'segmentation violation' "${LOG}gf_${r}.log" || [ ! -s "$fo" ]; then
    echo "[FAIL] $r  (see ${LOG}gf_${r}.log)"; rm -f "$fo"
  else
    grep -q "CATIMA material model: MSC ON, straggling ON" "${LOG}gf_${r}.log" \
      || echo "[WARN] $r: CATIMA line missing from the log -- backend may not be active"
    grep -q "dE/dx from CATIMA" "${LOG}gf_${r}.log" \
      || echo "[WARN] $r: CATIMA dE/dx banner missing -- fell back to the table"
    touch "${fo}.done"; echo "[ok] $r  $(date '+%H:%M:%S')"
  fi
}
export -f one

echo "=== CATIMA (d,p) production: $(echo $RUNS | wc -w) runs, $NPAR parallel -> $GF ==="
echo "    gate $GATE   table $TAB"
printf '%s\n' $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "=== done: $(ls "$GF"/*_genfitter_p.root 2>/dev/null | wc -l) fit files in $GF ==="
