#!/usr/bin/env bash
# Build the inelastic decks for the CONSISTENT optical-model test.
#
# The point of the test: the luminosity is measured on the elastic against a calculated
# sigma_el, so it carries the optical model; but the SAME potential also generates the inelastic
# distorted waves. Changing potential therefore moves L and the DWBA cross section in the same
# direction and the two partly cancel in the extracted B(EL). Quoting the L spread alone
# overstates the systematic. To measure the cancellation each potential has to be used in BOTH
# places, which is what these decks are for.
#
# The potential parameters are not typed in by hand -- they are taken from Ptolemy's own library
# via `--create-infile`, which writes the fully expanded deck. Reading them back off the printed
# potential table instead loses the imaginary spin-orbit and the Coulomb radius (KD03 has
# vsoi = -0.047 and rc0 = 1.478, both absent from that table).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PTOLEMY="${PTOLEMY:-/home/yassid/PtolemyCpp/ptolemy}"
mkdir -p "$HERE/expanded"
# level tag -> "Jpi Ex"
declare -A LEVELS=( [6094_1m]="1- 6.094" [6728_3m]="3- 6.728" [7012_2p]="2+ 7.012" )
for P in K V G P M; do
  printf '14C(p,p)14C  0+  none  0+  0.000  11.581MeV  %s\n' "$P" > "$HERE/expanded/el_$P.dwba"
  "$PTOLEMY" --create-infile "$HERE/expanded/el_$P.in" "$HERE/expanded/el_$P.dwba" >/dev/null 2>&1
  # the potential block is everything between the reference comment and ELASTIC SCATTERING
  POT=$(sed -n '/^\$.*http\|^\$.*Bardayan/,/^ELASTIC SCATTERING/p' "$HERE/expanded/el_$P.in" \
        | grep -E '^(v|vi|vsi|vso|vsoi) ')
  [ -n "$POT" ] || { echo "could not extract the potential for $P"; exit 1; }
  for L in "${!LEVELS[@]}"; do
    read -r JPI EX <<< "${LEVELS[$L]}"
    cat > "$HERE/inputs/omp_${P}_${L}.in" <<EOF
HEADER: 14C(P,P')14C* $JPI $EX MeV, OMP=$P, INELOCA1  [consistent OMP test]
JBIGA=0+
REACTION: 14C(P, P)14C($JPI $EX) ELAB=11.581
BELX = 0.1
PARAMETERSET INELOCA1
LMAX=30
INCOMING
r0target
$POT
;
OUTGOING
;
ANGLEMAX=180 ANGLESTEP=1
;
RETURN
EOF
  done
done
echo "  built $(ls "$HERE"/inputs/omp_*.in | wc -l) decks (5 potentials x 3 levels)"
