#!/usr/bin/env bash
# Browser explorer for 16C(d,t)15C reconstructed from the ARC-FIT Brho, not the genfit momentum.
#
# WHY: genfit is unusable for this channel. A 1 MeV triton has a 144 mm range and an 81 mm
# cyclotron radius in 2.85 T, so it stops before completing an arc and the fit reads the
# barely-curved fragment as high momentum -- 2.3x the energy the rigidity implies, rising to
# 25x for the shortest tracks. Ex carries a x5 lever arm on KE at forward angles, so that
# alone destroys the spectrum. The circle fit to the first arc has no such failure mode.
#
# The cache therefore holds  ke = KE(Brho, triton)  and  theta = Spyral polar (which IS
# theta_lab in Spyral's convention; ATTPCROOT's AtSpyralPID reports the supplement).
#
# The `chi2ndf` slider is REPURPOSED: there is no fit, so the column holds 100/sqrt(dE/dx).
# Lower = more ionising = more triton-like. sqrt(dE/dx) > 13  <=>  chi2ndf < 7.7.
set -o pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EBEAM="${1:-184.17}"
CACHE=/mnt/f/a1975/caches/dt_kin_spyral_brho.root
TMP="${DTTMP:-/mnt/f/a1975/caches/.explorer}"   # stable; the old value was a dead session scratchpad
mkdir -p "$TMP"
OUT="$HOME/a1975_C16_dt_brho_explorer.html"
MK="$HERE/../../../a2091/UKF/pp/make_explorer_html.C"
[[ -f "$CACHE" ]] || { echo "ERROR: no cache at $CACHE"; exit 1; }
set +u; source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh" >/dev/null 2>&1; set -u
root -b -l -q "$HERE/mkexp_pp.C(\"$CACHE\",\"$TMP/exp_dt_brho.root\",1e9,0,1e9)"
root -b -l -q "$MK(\"$TMP/exp_dt_brho.root\",\"$OUT\",\"16C(d,t)15C  from B#rho (no genfit)\",$EBEAM,16.0147013,2.0135532,3.01550072,15.0105993,16,\"0:g.s.,0.740:5/2+,1.218:Sn,3.103,4.220,4.657\",\"\")"
[[ -s "$OUT" ]] || { echo "ERROR: explorer not written"; exit 1; }
python3 "$HERE/add_keoff.py" "$OUT"
cp "$OUT" /mnt/c/Users/Yassid/Desktop/ 2>/dev/null
echo "explorer -> $OUT"
