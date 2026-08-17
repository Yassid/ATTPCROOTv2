#!/usr/bin/env bash
# Browser explorer for 16C(d,t)15C carrying the CATIMA material-effects production alongside
# the Spyral solver and the production in use today, so the three compare on identical axes:
#
#     CATIMA          genfit + CATIMA matFX, after the TGeo navigator fix
#                                                     /mnt/f/a1975/caches/dt_kin_catima.root
#     Spyral_solver   InterpSolverPhase, RK trajectory fit
#                                                     /mnt/f/a1975/caches/dt_kin_spyral.root
#     genfit_matOFF   genfit, material effects off    /mnt/f/a1975/caches/dt_kin_dv1104.root
#
# CATIMA is first so it is the page's default set.
#
# COLUMN CHOICE, not uniform across the sets and not to be "tidied":
#   CATIMA  ke/theta -- GetKinematicsXtr, the back-extrapolated (corrected) slot
#   Spyral  ke/theta -- the solver's own result, already at the vertex
#   genfit  ke/theta -- GetKinematicsXtr as well
# (Only the UKF differs, keeping its corrected value in kefit/thetafit. See
#  open_explorer_dt_spyral.sh, which carries the UKF instead of CATIMA.)
#
# IC GATE 900-1400 ON THE TWO ATTPCROOT SETS, and it is required. Spyral is already IC-gated
# because InterpSolverPhase takes ic_min_val/ic_max_val and only solves what passes, while the
# ATTPCROOT caches merely STORE ic. Leaving them ungated would put a gated set beside two still
# carrying the ~38% beam contaminant.
#
# TWO DIFFERENCES FROM open_explorer_dt_spyral.sh, both deliberate:
#
#   chi2max = 1e8 rather than the conventional 1e9. mkexp_pp now drops collapsed fits
#   (chi2ndf >= 1e9, written when ndf <= 0) unconditionally, so this is belt-and-braces and
#   every sibling script is safe too -- but 1e8 also keeps the intent visible at the call site.
#
#   TMP is a stable directory, not a session scratchpad. The older open_explorer_* scripts point
#   at /tmp/claude-1000/.../b025789e-.../scratchpad, which no longer exists, so they fail today.
#
#   ./open_explorer_dt_catima.sh [Ebeam]
set -o pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EBEAM="${1:-184.17}"
TMP="/mnt/f/a1975/caches/.explorer"
OUT="$HOME/a1975_C16_dt_catima_explorer.html"
MK="$HERE/../../../a2091/UKF/pp/make_explorer_html.C"

C_CAT="/mnt/f/a1975/caches/dt_kin_catima.root"
C_SPY="/mnt/f/a1975/caches/dt_kin_spyral.root"
C_OFF="/mnt/f/a1975/caches/dt_kin_dv1104.root"
mkdir -p "$TMP"
set +u; source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh" >/dev/null 2>&1; set -u

for c in "$C_CAT" "$C_SPY" "$C_OFF"; do
  [[ -s "$c" ]] || { echo "ERROR: missing cache $c"; exit 1; }
done
[[ -s "$MK" ]] || { echo "ERROR: missing $MK"; exit 1; }

root -b -l -q "$HERE/mkexp_pp.C(\"$C_CAT\",\"$TMP/exp_dt_cat.root\",1e8,900,1400,\"ke\",\"theta\")"
root -b -l -q "$HERE/mkexp_pp.C(\"$C_SPY\",\"$TMP/exp_dt_spy.root\",1e8,0,1e9,\"ke\",\"theta\")"
root -b -l -q "$HERE/mkexp_pp.C(\"$C_OFF\",\"$TMP/exp_dt_off.root\",1e8,900,1400,\"ke\",\"theta\")"

# 15C: g.s., the 0.740 5/2+ (the only bound excited state), Sn at 1.218, then the unbound structures
root -b -l -q "$MK(\"$TMP/exp_dt_cat.root\",\"$OUT\",\"16C(d,t)15C  CATIMA vs Spyral vs genfit matFX-off\",$EBEAM,16.0147013,2.0135532,3.01550072,15.0105993,16,\"0:g.s.,0.740:5/2+,1.218:Sn,3.103,4.220,4.657\",\"$TMP/exp_dt_spy.root\",\"$TMP/exp_dt_off.root\",\"CATIMA,Spyral_solver,genfit_matOFF\")"
[[ -s "$OUT" ]] || { echo "ERROR: explorer not written"; exit 1; }
python3 "$HERE/add_keoff.py" "$OUT"
cp "$OUT" /mnt/c/Users/Yassid/Desktop/ 2>/dev/null && echo "copied to the Windows Desktop"
echo "explorer -> $OUT"
