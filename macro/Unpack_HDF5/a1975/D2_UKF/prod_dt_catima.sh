#!/usr/bin/env bash
# a1975 16C(d,t)15C genfit production at dv 1.10424, WITH material effects, using the CATIMA
# backend. Supersedes the genfit arm of fit_dt_dv1104.sh, which runs matEffects = kFALSE.
#
# WHAT CHANGED FROM fit_dt_dv1104.sh, AND WHY EACH ONE
#
#   matEffects kFALSE -> kTRUE.  The point of the exercise. It was off because genfit's own
#       material models leave 62% of triton fits with ndf < 0 -- collapsed fits that still carry
#       kinematics and so are invisible unless ndf is checked. With CATIMA that is 3%.
#
#   catimaMSC = catimaStraggling = kTRUE.  Needs a GenFit built -DGENFIT_USE_CATIMA=ON. Both are
#       required: the collapse is a CONJUNCTION, and removing either one (or the table below)
#       puts it back to ~50%. Straggling is the term that actually drives it; the MSC model on
#       its own moves 62% -> 56%.
#
#   eLossTable = triton_D2_300torr.txt.  Below beta*gamma = 0.05 (KE = 3.5 MeV for a triton)
#       genfit applies NO stopping power without a curve loaded, and the (d,t) low branch is
#       0.8-6 MeV. Regenerated 2026-08-16 with deuterium's real mass 2.014 rather than a round 2,
#       which had been a flat +0.70% on every point.
#
#   manualElossDensity 6.61e-5 -> 0.  THE ONE THAT IS EASY TO GET WRONG. With material effects
#       ON, extrapolateToLine already transports the state back through the gas with genfit's own
#       RK4 dE/dx, so the vertex-gap loss is ALREADY in KE_xtr. SetManualELoss then adds
#       GetEnergyLoss(KE, vtxGapCm) on top of it -- it is documented as being for the
#       material-effects-OFF case, where the extrapolation is geometric and leaves |p| untouched,
#       but nothing in the code enforces that. Passing 0 disables ONLY SetManualELoss:
#       SetELossHybrid falls back to 6.61e-5 internally and the range constraint is separate.
#       Measured cost of leaving it in: +5.25% on the low branch of KE_xtr, which is the slot the
#       Ex analysis reads, and a 14% wider low-branch distribution.
#
#   matFallback kTRUE -> kFALSE.  A throwing track is NOT silently refitted with
#       setNoEffects(true) and kept. This keeps the production a single physics model and
#       directly comparable to the run_0031 arms. The alternative (kTRUE) would preserve a few
#       percent more yield but every consumer would then have to filter on
#       GetMatEffectsFallback(), and not filtering is exactly how a blended spectrum happens.
#
# UNCHANGED: gate, geometry, par file, beam energy, gas density, iteration counts, theta window,
# measSigma. Same reco input as the dv 1.10424 production.
#
#   ./prod_dt_catima.sh [nparallel] [runs...]
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
export GF=/mnt/f/a1975/gf_dt_catima/
export LOG=/mnt/f/a1975/logs_catima/
export PAR=ATTPC.a1975_deuterium_dv1104.par
export GATE=pid/triton_d2_dv1104.json
export TAB=triton_D2_300torr.txt
mkdir -p "$GF" "$LOG"

one() {
  n="$1"; r="run_${n}"
  rc="${REC}${r}_multifit_reco.root"
  [ -s "$rc" ] || { echo "[noreco] $r"; return 0; }
  fo="${GF}${r}_multifit_genfitter_t.root"
  [ -s "${fo}.done" ] && { echo "[have] $r"; return 0; }
  root -l -b -q "fitGenfitter_a1975_deuterium.C(\"${r}_multifit\",-1,\"$REC\",\"\",\"$GF\",\
-2.85,2,5,\"$GATE\",4.0,10.0,170.0,kTRUE,kTRUE,1000010030,3.01550072,1,\"t\",\"_reco\",\
\"ATTPC_D300torr_v2_geomanager.root\",kTRUE,0,2,\"$PAR\",kFALSE,kFALSE,kFALSE,\"$TAB\",kTRUE,kTRUE)" \
    > "${LOG}gf_${r}.log" 2>&1
  if grep -qi 'segmentation violation' "${LOG}gf_${r}.log" || [ ! -s "$fo" ]; then
    echo "[FAIL] $r  (see ${LOG}gf_${r}.log)"; rm -f "$fo"
  else
    # a run whose CATIMA backend silently did not engage is worth catching here, not at analysis
    grep -q "CATIMA material model: MSC ON, straggling ON" "${LOG}gf_${r}.log" \
      || echo "[WARN] $r: CATIMA line missing from the log -- backend may not be active"
    touch "${fo}.done"; echo "[ok] $r  $(date '+%H:%M:%S')"
  fi
}
export -f one

echo "=== CATIMA (d,t) production: $(echo $RUNS | wc -w) runs, $NPAR parallel -> $GF ==="
printf '%s\n' $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "=== done: $(ls "$GF"/*_genfitter_t.root 2>/dev/null | wc -l) fit files in $GF ==="
