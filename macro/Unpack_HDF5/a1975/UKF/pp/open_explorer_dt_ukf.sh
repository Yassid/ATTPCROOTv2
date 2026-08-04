#!/usr/bin/env bash
# Browser explorer for a1975 16C(d,t)15C: UKF against GENFIT on the SAME dv=1.136 point clouds.
#
# This is the page's fitter switch used for what it was built for. Both slots are the same 8
# runs, same PID-gated triton chain, same reconstruction -- the only difference is the fitter,
# and with it the energy loss: AtFitterUKF back-extrapolates to the beam axis (AtFitterUKF.cxx:713)
# and propagates with an AtELossModel, while the genfit path here runs matEffects = kFALSE with
# no dE/dx table, so it reports the state at the first measurement point with no energy restored.
#
# Measured on these two caches (g.s. / 3.103 / spacing, truth 0 / 3.103 / 3.103):
#   genfit  -0.040 / +2.879 / 2.919   sigma 0.252
#   UKF     +0.019 / +3.211 / 3.192   sigma 0.155
# The UKF gas density is 1.322e-4 g/cm3 (D2, 600 torr, 293 K), NOT the macro default 9.0e-5 --
# at the default the spacing collapses to 2.805, worse than genfit.
#
#   ./open_explorer_dt_ukf.sh          # Ebeam 180
set -o pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EBEAM="${1:-180}"
SCR="${UKFSCR:-/tmp/claude-1000/-home-yassid/9b2b9ba6-e9b5-4d13-af2f-ad9750f083fc/scratchpad}"
CUKF="${SCR}/ukf/dt_ukf_d132.root"; CGF="${SCR}/dv1136/dt_1136.root"
TMP="${UKFTMP:-/tmp}"; OUT="$HOME/a1975_C16_dt_ukf_vs_genfit.html"
MK="$HERE/../../../a2091/UKF/pp/make_explorer_html.C"
for f in "$MK" "$CUKF" "$CGF"; do [[ -f "$f" ]] || { echo "ERROR: missing $f"; exit 1; }; done

set +u; source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh" >/dev/null 2>&1; set -u

# D2 beam gate [900,1300], not mkexp_pp's H2 default [950,1350]
root -b -l -q "$HERE/mkexp_pp.C(\"$CUKF\",\"$TMP/exp_ukf.root\",1e9,900,1300)"
root -b -l -q "$HERE/mkexp_pp.C(\"$CGF\",\"$TMP/exp_gf.root\",1e9,900,1300)"
# first slot is labelled UKF and second GENFIT by the builder, which is correct here for once,
# so no relabelling step -- unlike the dv-scan page where the slots hold drift velocities
root -b -l -q "$MK(\"$TMP/exp_ukf.root\",\"$OUT\",\"16C(d,t)15C  UKF vs GENFIT\",$EBEAM,16.0147013,2.0135532,3.01550072,15.0105993,16,\"0:g.s.,0.740:5/2+,1.218:Sn,3.103,4.220,4.657\",\"$TMP/exp_gf.root\")"
[[ -s "$OUT" ]] || { echo "ERROR: explorer not written"; exit 1; }

if grep -qi microsoft /proc/version 2>/dev/null; then
  WINHOME="/mnt/c/Users/$(ls /mnt/c/Users | grep -viE 'public|default|all users' | head -1)"
  BROWSER="/mnt/c/Program Files/Google/Chrome/Application/chrome.exe"
  [[ -x "$BROWSER" ]] || BROWSER="/mnt/c/Program Files (x86)/Microsoft/Edge/Application/msedge.exe"
  cp "$OUT" "$WINHOME/$(basename "$OUT")"
  WINPATH="$(wslpath -w "$WINHOME/$(basename "$OUT")" | sed 's|\\|/|g')"
  echo "opening file:///$WINPATH"
  nohup "$BROWSER" "file:///$WINPATH" >/dev/null 2>&1 &
fi
echo "wrote $OUT"
