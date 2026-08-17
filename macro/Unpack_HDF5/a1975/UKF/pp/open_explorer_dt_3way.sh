#!/usr/bin/env bash
# Browser explorer for 16C(d,t)15C at dv 1.10424, Ebeam 184.17 MeV, carrying THREE fitters
# side by side on one page with a selector button each:
#
#     genfit_matOFF   genfit, material effects off      /mnt/f/a1975/caches/dt_kin_dv1104.root
#     genfit_matON    genfit, material effects on       /mnt/f/a1975/caches/dt_kin_maton.root
#     UKF             AtFitterUKF (RK4 + dE/dx)         /mnt/f/a1975/caches/dt_kin_ukf.root
#
# WHY THE UKF SET USES DIFFERENT COLUMNS.  The two fitters swap their kinematics slots:
# for genfit the corrected (back-extrapolated) energy is GetKinematicsXtr, which ex_dt_a1975.C
# writes as `ke`/`theta`; for the UKF the corrected (at-vertex) energy is GetKinematics, which
# ex_dt_a1975.C writes as `kefit`/`thetafit`. Feeding the page `ke` for both would silently
# compare genfit's corrected energy against the UKF's RAW first-cluster energy. Hence the
# keCol/thCol arguments to mkexp_pp below -- do not "tidy" them into one call.
#
# The matON production is fitted with the matFX fallback OFF, so that set is pure material
# effects rather than a blend of matFX and no-matFX retries.
#
# NO CUTS ARE PRE-APPLIED: chi2 -> 1e9, IC window 0..1e9. Every selection is made in the page.
#
#   ./open_explorer_dt_3way.sh [Ebeam]
set -o pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EBEAM="${1:-184.17}"
TMP="${DTTMP:-/mnt/f/a1975/caches/.explorer}"   # stable; the old value was a dead session scratchpad
mkdir -p "$TMP"
OUT="$HOME/a1975_C16_dt_3way_explorer.html"
MK="$HERE/../../../a2091/UKF/pp/make_explorer_html.C"

C_OFF="/mnt/f/a1975/caches/dt_kin_dv1104.root"
C_ON="/mnt/f/a1975/caches/dt_kin_maton.root"
C_UKF="/mnt/f/a1975/caches/dt_kin_ukf.root"
set +u; source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh" >/dev/null 2>&1; set -u

for c in "$C_OFF" "$C_ON" "$C_UKF"; do
  [[ -s "$c" ]] || { echo "ERROR: missing cache $c"; exit 1; }
done

root -b -l -q "$HERE/mkexp_pp.C(\"$C_OFF\",\"$TMP/exp_dt_off.root\",1e9,0,1e9,\"ke\",\"theta\")"
root -b -l -q "$HERE/mkexp_pp.C(\"$C_ON\",\"$TMP/exp_dt_on.root\",1e9,0,1e9,\"ke\",\"theta\")"
# UKF: corrected energy lives in the OTHER slot -- see the note above
root -b -l -q "$HERE/mkexp_pp.C(\"$C_UKF\",\"$TMP/exp_dt_ukf.root\",1e9,0,1e9,\"kefit\",\"thetafit\")"

# 15C levels: g.s., 0.740 (5/2+), Sn at 1.218, then the unbound structures.
root -b -l -q "$MK(\"$TMP/exp_dt_off.root\",\"$OUT\",\"16C(d,t)15C  dv 1.10424  --  3 fitters\",$EBEAM,16.0147013,2.0135532,3.01550072,15.0105993,16,\"0:g.s.,0.740:5/2+,1.218:Sn,3.103,4.220,4.657\",\"$TMP/exp_dt_on.root\",\"$TMP/exp_dt_ukf.root\",\"genfit_matOFF,genfit_matON,UKF\")"
[[ -s "$OUT" ]] || { echo "ERROR: explorer not written"; exit 1; }
python3 "$HERE/add_keoff.py" "$OUT"
cp "$OUT" /mnt/c/Users/Yassid/Desktop/ 2>/dev/null && echo "copied to the Windows Desktop"
echo "explorer -> $OUT"
