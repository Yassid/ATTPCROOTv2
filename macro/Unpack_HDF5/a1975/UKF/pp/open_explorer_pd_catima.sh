#!/usr/bin/env bash
# Browser explorer for 16C(p,d)15C carrying the CATIMA material-effects production beside the
# production the adopted result was built on, so the two compare on identical axes:
#
#     CATIMA          genfit + CATIMA matFX, after the TGeo navigator fix
#                                              /mnt/f/a1975/caches/pd_kin_catima.root   (30768)
#     genfit_matOFF   genfit, material effects off
#                                              /mnt/f/a1975_C16_pd_results/pd_kin.root  (28346)
#
# CATIMA is first so it is the page's default set.
#
# THE BASELINE IS pd_kin.root, NOT UKF/deuteron_kin.root. The latter is a June 10 file with 19473
# tracks and a different selection; pairing it with the CATIMA cache invents a +58% statistics
# gain that is really a selection difference. Against pd_kin.root the honest figure is +8.5%.
#
# COLUMN CHOICE. Both sets use ke/theta because that is the ONLY kinematics these caches carry:
# cache_pd_run.C:45 reads GetKinematics(), which for genfit is the RAW slot -- the fit at the
# first measurement point, NOT back-extrapolated to the vertex. The (d,t) caches carry both slots
# and their explorers use the corrected one, so DO NOT read absolute Ex here against a (d,t) page.
# The comparison on this page is unaffected: both sets are built the same way, so the only thing
# that differs between them is the fitter.
#
# COLLAPSED FITS. This is the first (p,d) production with material effects ON, which is precisely
# when the chi2ndf = 1e9 sentinel (written when ndf <= 0) stops being negligible -- on (d,t) with
# Highland it was 60.4% of everything the page displayed. mkexp_pp drops it unconditionally and
# reports the count, so watch the "collapsed" figure in the output of the two calls below.
#
# IC GATE 950-1350, mkexp_pp's own default, which is how the existing (p,d) page was built. Passed
# explicitly here so the number is visible rather than inherited.
#
# EBEAM defaults to 185, the value in the saved explorer config behind the adopted (p,d) cross
# sections. The older open_explorer_pd.sh defaults to 195.5; that disagreement is not resolved
# here, so pass the value you want as the first argument.
#
# NOTE: no `set -e` -- this ROOT build segfaults in TROOT::EndOfProcessCleanups AFTER writing its
# output, so every root call returns non-zero on success.
#
#   ./open_explorer_pd_catima.sh [Ebeam]
set -o pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EBEAM="${1:-185}"
TMP="${PDTMP:-/mnt/f/a1975/caches/.explorer}"
OUT="$HOME/a1975_C16_pd_catima_explorer.html"
MK="$HERE/../../../a2091/UKF/pp/make_explorer_html.C"

C_CAT="/mnt/f/a1975/caches/pd_kin_catima.root"
C_OFF="/mnt/f/a1975_C16_pd_results/pd_kin.root"
mkdir -p "$TMP"
set +u; source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh" >/dev/null 2>&1; set -u

for c in "$C_CAT" "$C_OFF"; do
  [[ -s "$c" ]] || { echo "ERROR: missing cache $c"; exit 1; }
done
[[ -s "$MK" ]] || { echo "ERROR: missing $MK"; exit 1; }

root -b -l -q "$HERE/mkexp_pp.C(\"$C_CAT\",\"$TMP/exp_pd_cat.root\",1e8,950,1350,\"ke\",\"theta\")"
root -b -l -q "$HERE/mkexp_pp.C(\"$C_OFF\",\"$TMP/exp_pd_off.root\",1e8,950,1350,\"ke\",\"theta\")"

# 16C(p,d)15C: beam 16C, target p, ejectile d, residual 15C, beamA=16.
# 15C levels: g.s., the 0.740 5/2+ (the only bound excited state), Sn at 1.218, then unbound.
root -b -l -q "$MK(\"$TMP/exp_pd_cat.root\",\"$OUT\",\"16C(p,d)15C  CATIMA vs genfit matFX-off\",$EBEAM,16.0147013,1.00782503,2.0135532,15.0105993,16,\"0:g.s.,0.740:5/2+,1.218:Sn,3.103,4.220,4.657\",\"$TMP/exp_pd_off.root\",\"\",\"CATIMA,genfit_matOFF\")"
[[ -s "$OUT" ]] || { echo "ERROR: explorer not written"; exit 1; }
python3 "$HERE/add_keoff.py" "$OUT"

# Claude Code has no rendering surface, so the page has to be handed to a Windows browser.
if grep -qi microsoft /proc/version 2>/dev/null; then
  cp "$OUT" /mnt/c/Users/Yassid/Desktop/ 2>/dev/null && echo "copied to the Windows Desktop"
  WINPATH="$(wslpath -w "/mnt/c/Users/Yassid/Desktop/$(basename "$OUT")" 2>/dev/null)"
  [[ -n "$WINPATH" ]] && cmd.exe /c start "" "$WINPATH" >/dev/null 2>&1 && echo "opened in the default browser"
fi
echo "explorer -> $OUT"
