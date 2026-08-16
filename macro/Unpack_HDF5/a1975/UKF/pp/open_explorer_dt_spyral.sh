#!/usr/bin/env bash
# Browser explorer for 16C(d,t)15C carrying the SPYRAL SOLVER result alongside the two
# ATTPCROOT fitters, so the three can be compared on identical axes in one page:
#
#     Spyral_solver   InterpSolverPhase, RK trajectory fit   /mnt/f/a1975/caches/dt_kin_spyral.root
#     genfit_matOFF   genfit, material effects off           /mnt/f/a1975/caches/dt_kin_dv1104.root
#     UKF             AtFitterUKF                            /mnt/f/a1975/caches/dt_kin_ukf.root
#
# Spyral is first so it is the page's default set.
#
# COLUMN CHOICE, which is not uniform across the three and must not be "tidied":
#   Spyral  ke/theta   -- the solver's own result, already at the vertex
#   genfit  ke/theta   -- GetKinematicsXtr, the back-extrapolated (corrected) slot
#   UKF     kefit/thetafit -- GetKinematics, which is the CORRECTED slot for the UKF; the two
#           ATTPCROOT fitters store their corrected value in OPPOSITE slots.
#
# The Spyral cache is built by mk_spyral_cache.C from a text dump rather than by uproot:
# uproot has written an RNTuple here before, which ROOT 6.26 cannot read.
#
# IC GATE 900-1400 IS APPLIED TO ALL THREE SETS, and it has to be. The Spyral set is already
# IC-gated because InterpSolverPhase takes ic_min_val/ic_max_val and only solves what passes,
# whereas the ATTPCROOT caches STORE ic and never cut on it (cache_dt_dv1104.sh keeps the beam
# gate tunable at analysis time). Leaving the ATTPCROOT sets ungated would put a gated set next
# to two that still carry the ~38% beam contaminant, and the selector would be comparing
# different beams. chi2 stays open at 1e9; that selection is made in the page.
#
#   ./open_explorer_dt_spyral.sh [Ebeam]
set -o pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EBEAM="${1:-184.17}"
TMP="/tmp/claude-1000/-home-yassid/b025789e-3ded-4e0d-8d25-35b205d047eb/scratchpad"
OUT="$HOME/a1975_C16_dt_spyral_explorer.html"
MK="$HERE/../../../a2091/UKF/pp/make_explorer_html.C"

C_SPY="/mnt/f/a1975/caches/dt_kin_spyral.root"
C_OFF="/mnt/f/a1975/caches/dt_kin_dv1104.root"
C_UKF="/mnt/f/a1975/caches/dt_kin_ukf.root"
set +u; source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh" >/dev/null 2>&1; set -u

for c in "$C_SPY" "$C_OFF" "$C_UKF"; do
  [[ -s "$c" ]] || { echo "ERROR: missing cache $c"; exit 1; }
done

root -b -l -q "$HERE/mkexp_pp.C(\"$C_SPY\",\"$TMP/exp_dt_spy.root\",1e9,0,1e9,\"ke\",\"theta\")"
root -b -l -q "$HERE/mkexp_pp.C(\"$C_OFF\",\"$TMP/exp_dt_off.root\",1e9,900,1400,\"ke\",\"theta\")"
root -b -l -q "$HERE/mkexp_pp.C(\"$C_UKF\",\"$TMP/exp_dt_ukf.root\",1e9,900,1400,\"kefit\",\"thetafit\")"

# 15C: g.s., the 0.740 5/2+ (the only bound excited state), Sn at 1.218, then the unbound structures
root -b -l -q "$MK(\"$TMP/exp_dt_spy.root\",\"$OUT\",\"16C(d,t)15C  Spyral solver vs ATTPCROOT\",$EBEAM,16.0147013,2.0135532,3.01550072,15.0105993,16,\"0:g.s.,0.740:5/2+,1.218:Sn,3.103,4.220,4.657\",\"$TMP/exp_dt_off.root\",\"$TMP/exp_dt_ukf.root\",\"Spyral_solver,genfit_matOFF,UKF\")"
[[ -s "$OUT" ]] || { echo "ERROR: explorer not written"; exit 1; }
python3 "$HERE/add_keoff.py" "$OUT"
cp "$OUT" /mnt/c/Users/Yassid/Desktop/ 2>/dev/null && echo "copied to the Windows Desktop"
echo "explorer -> $OUT"
