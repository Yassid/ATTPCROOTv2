#!/usr/bin/env bash
# Browser explorer for 16C(d,t)15C: does NATIVE CATIMA dE/dx change the spectrum, or only the
# number of fits that survive? Three sets on identical axes:
#
#     CATIMA_native   genfit, dE/dx from catima::dedx per step, the STEP's own material
#                                              /mnt/f/a1975/caches/dt_kin_cateloss.root  (36119 usable)
#     CATIMA_table    genfit, dE/dx from the frozen ASCII curve x one global density
#                                              /mnt/f/a1975/caches/dt_kin_catima.root    (35601 usable)
#     Spyral_solver   InterpSolverPhase, RK trajectory fit -- the ADOPTED (d,t) energies
#                                              /mnt/f/a1975/caches/dt_kin_spyral.root
#
# The first two differ in ONE flag (catimaELoss) and nothing else, so anything you see between
# them is the dE/dx source. Spyral is the reference: no genfit/UKF energy resolves the two low
# 15C levels, so the solver remains the source of the adopted numbers regardless of this test.
#
# WHAT THE FITS ALREADY SAY, so the page can be judged against a prediction rather than a hope:
# collapse 2.51% -> 1.08% (+518 usable tracks), 3.103 peak 2.872 -> 2.866 (-0.006 +- 0.018) and
# sigma 0.247 -> 0.252 (+0.005 +- 0.021). i.e. MORE STATISTICS, SAME SPECTRUM. If the two curves
# look different in shape or position here, something is wrong with the page, not with the physics.
# Both remain ~0.24 MeV below the known 3.103 -- that deficit is NOT the energy loss.
#
# COLUMN CHOICE, not uniform and not to be "tidied": ke/theta everywhere, but for the genfit sets
# that is GetKinematicsXtr (back-extrapolated) while for Spyral it is the solver's own vertex
# result. See reference_genfit_ukf_slots.
#
# IC GATE 900-1400 on the two genfit sets and required: Spyral is already gated at solve time
# (InterpSolverPhase takes ic_min_val/ic_max_val) while the ATTPCROOT caches only STORE ic, so
# leaving them ungated would put a gated set beside two carrying the ~38% beam contaminant.
#
# chi2max = 1e8, not the conventional 1e9: 1e9 is ALSO the collapsed-fit sentinel that
# ex_dt_a1975 writes when ndf <= 0, so `chi2ndf > 1e9` is false and collapsed fits sail through.
# mkexp_pp drops them unconditionally anyway; 1e8 keeps the intent visible at the call site.
# That matters more here than usual -- the whole point of this comparison is a collapse rate.
#
#   ./open_explorer_dt_eloss.sh [Ebeam]
set -o pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EBEAM="${1:-184.17}"
TMP="${DTTMP:-/mnt/f/a1975/caches/.explorer}"
OUT="$HOME/a1975_C16_dt_bethebloch_explorer.html"
MK="$HERE/../../../a2091/UKF/pp/make_explorer_html.C"

C_NEW="/mnt/f/a1975/caches/dt_kin_catfull.root"
C_OLD="/mnt/f/a1975/caches/dt_kin_cateloss.root"
C_SPY="/mnt/f/a1975/caches/dt_kin_catima.root"
mkdir -p "$TMP"
set +u; source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh" >/dev/null 2>&1; set -u

for c in "$C_NEW" "$C_OLD" "$C_SPY"; do
  [[ -s "$c" ]] || { echo "ERROR: missing cache $c"; exit 1; }
done
[[ -s "$MK" ]] || { echo "ERROR: missing $MK"; exit 1; }

root -b -l -q "$HERE/mkexp_pp.C(\"$C_NEW\",\"$TMP/exp_dt_full.root\",1e8,900,1400,\"ke\",\"theta\")"
root -b -l -q "$HERE/mkexp_pp.C(\"$C_OLD\",\"$TMP/exp_dt_low.root\",1e8,900,1400,\"ke\",\"theta\")"
root -b -l -q "$HERE/mkexp_pp.C(\"$C_SPY\",\"$TMP/exp_dt_tab.root\",1e8,900,1400,\"ke\",\"theta\")"

# 15C: g.s., the 0.740 5/2+ (the only bound excited state), Sn at 1.218, then the unbound structures
root -b -l -q "$MK(\"$TMP/exp_dt_full.root\",\"$OUT\",\"16C(d,t)15C  native CATIMA dE/dx vs table vs Spyral\",$EBEAM,16.0147013,2.0135532,3.01550072,15.0105993,16,\"0:g.s.,0.740:5/2+,1.218:Sn,3.103,4.220,4.657\",\"$TMP/exp_dt_low.root\",\"$TMP/exp_dt_tab.root\",\"CATIMA_full,CATIMA_low_BB,table_BB\")"
[[ -s "$OUT" ]] || { echo "ERROR: explorer not written"; exit 1; }
python3 "$HERE/add_keoff.py" "$OUT"

# Claude Code has no rendering surface, so hand the page to a Windows browser.
if grep -qi microsoft /proc/version 2>/dev/null; then
  cp "$OUT" /mnt/c/Users/Yassid/Desktop/ 2>/dev/null && echo "copied to the Windows Desktop"
  WINPATH="$(wslpath -w "/mnt/c/Users/Yassid/Desktop/$(basename "$OUT")" 2>/dev/null)"
  [[ -n "$WINPATH" ]] && cmd.exe /c start "" "$WINPATH" >/dev/null 2>&1 && echo "opened in the default browser"
fi
echo "explorer -> $OUT"
