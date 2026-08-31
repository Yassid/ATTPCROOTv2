#!/usr/bin/env bash
# Build one 46Ar(3He,d) parameter file per magnetic field from the Magboltz 300 torr scan.
#
#   ./make_3Hed_pars.sh
#
# WHY THIS EXISTS. Until 2026-08-30 the 2.85 T and 3.8 T pars were byte-identical except for the
# BField line: DriftVelocity 1.30 cm/us and CoefL = CoefT = 0.0009 cm^2/us were PLACEHOLDERS
# carried over from a1975 H2 and held FIXED across fields. A field comparison run that way varies
# only the magnetic curvature, while transverse diffusion -- which really does fall with B through
# the omega*tau suppression -- stays frozen. That biases the comparison AGAINST the higher field.
#
# UNITS. Magboltz reports DL/DT in um/sqrt(cm); AtClusterize wants cm^2/us, because it forms
# sigma = sqrt(2 * Coef * driftTime). Equating sigma^2 = D^2 * z with 2 * Coef * z / v_d gives
#
#       Coef [cm^2/us] = D[cm/sqrt(cm)]^2 * v_d[cm/us] / 2,     D[cm] = D[um] / 1e4
#
# Sanity check on the OLD placeholders: 0.0009 and v_d 1.30 invert to D_T = 372 um/sqrt(cm),
# a believable H2 number -- so the conversion is consistent with what the chain already did.
#
# The E point is 50 V/cm, which is the EField 5000 V/m the par itself specifies. If you change
# EField, re-run the scan; do not interpolate this table in your head.
set -eo pipefail
SCAN=${SCAN:-/home/yassid/attpc_dv_3He/out300}
PREFIX=${PREFIX:-he3co2_300torr_op140}   # the operating-point scan (E ~ 97 V/cm); "he3co2_300torr" is the old 40-60 V/cm one
PARDIR=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/parameters
BASE=$PARDIR/ATTPC.46Ar_3Hed_sim.par
EPOINT=${EPOINT:-140}      # V/cm -- v_d 1.32 cm/us at 300 torr (E/p 0.465)
EFIELD_VM=${EFIELD_VM:-14000} # the same field in V/m, written into the par
SAMPRATE=${SAMPRATE:-6}   # MHz; 512 tb / 6 MHz = 85.3 us covers the 75.8 us drift with 13% margin

# THE SAMPLING RATE IS NOT A FREE NUMBER. AtDigiPar::GetTBTime() (AtParameter/AtDigiPar.cxx:20) is
# a hardcoded switch over the real GET rates and returns -1 ns for anything else. A par written
# with SamplingRate 5 parses fine, runs to completion, and yields ZERO fitted tracks because the
# whole time->z mapping is built on TBTime = -1. That cost an 85-minute sample on 2026-08-31.
case " 3 6 12 25 50 100 " in
   *" $SAMPRATE "*) ;;
   *) echo "SamplingRate $SAMPRATE MHz is NOT a GET rate (3/6/12/25/50/100) -- AtDigiPar would return TBTime = -1 ns"; exit 2 ;;
esac

[ -f "$BASE" ] || { echo "missing base par $BASE"; exit 2; }

emit() { # B_tesla  suffix
  local B=$1 sfx=$2
  local f="$SCAN/${PREFIX}_B${B}.txt"
  [ -s "$f" ] || { echo "MISSING SCAN $f -- run run_300torr_bscan.sh"; exit 1; }
  # columns: E  E/p  E/N  vd  vderr  DL[um]  DT[um]  lorentz  alpha  eta
  local row
  row=$(awk -v e="$EPOINT" '!/^#/ && ($1+0)>e-0.5 && ($1+0)<e+0.5 {print; exit}' "$f")
  [ -n "$row" ] || { echo "no E=$EPOINT row in $f"; exit 1; }
  read -r E EP EN VD VDERR DL DT LOR ALPHA ETA <<< "$row"
  local coefL coefT
  coefL=$(awk -v d="$DL" -v v="$VD" 'BEGIN{printf "%.6g", (d/1e4)*(d/1e4)*v/2}')
  coefT=$(awk -v d="$DT" -v v="$VD" 'BEGIN{printf "%.6g", (d/1e4)*(d/1e4)*v/2}')
  local out="$PARDIR/ATTPC.46Ar_3Hed_sim_${sfx}.par"
  awk -v B="$B" -v vd="$VD" -v cl="$coefL" -v ct="$coefT" -v dl="$DL" -v dt="$DT" -v ve="$VDERR" -v ef="$EFIELD_VM" -v sr="$SAMPRATE" '
    /^BField:Double_t/        {printf "BField:Double_t          %12s   # Tesla\n", B; next}
    /^DriftVelocity:Double_t/ {printf "DriftVelocity:Double_t   %12.4f   # cm/us  MAGBOLTZ, 3He+5%%CO2 300 torr 50 V/cm, B = %s T (+-%.2f%%)\n", vd, B, ve; next}
    /^CoefL:Double_t/         {printf "CoefL:Double_t           %12s   # cm^2/us  from D_L = %.1f um/sqrt(cm)\n", cl, dl; next}
    /^CoefT:Double_t/         {printf "CoefT:Double_t           %12s   # cm^2/us  from D_T = %.1f um/sqrt(cm)\n", ct, dt; next}
    /^EField:Double_t/        {printf "EField:Double_t          %12s   # V/m -- the field this v_d was computed AT. The old par said 5000 V/m\n", ef; next}
    /^SamplingRate:Int_t/     {printf "SamplingRate:Int_t       %12s   # MHz; 512 tb / %s MHz = %.1f us, and the 100 cm drift takes %.1f us\n", sr, sr, 512/sr, 100/vd; next}
    {print}
  ' "$BASE" > "$out"
  # THE DRIFT MUST FIT THE READOUT WINDOW. At the placeholder v_d this was never tested, and at
  # the true v_d with the original 512 tb / 6 MHz it FAILS: 165 us against an 85 us window, so
  # everything past z ~ 52 cm would have been silently lost. Assert it, do not assume it.
  awk -v v="$VD" -v sr="$SAMPRATE" 'BEGIN{
     t=100.0/v; w=512.0/sr;
     if (t>w) { printf "  *** DRIFT OVERRUNS THE WINDOW: %.1f us of drift against %.1f us of readout ***\n", t, w; exit 1 }
     printf "  drift %.1f us fits the %.1f us window (%.0f%% margin), %.2f mm per time bucket\n", t, w, 100*(w-t)/t, v/sr*10
  }' || exit 1
  printf "%-8s B=%-5s v_d=%8.4f cm/us  D_L=%7.1f  D_T=%7.1f um/sqrt(cm)  ->  CoefL=%-10s CoefT=%-10s  %s\n" \
         "$sfx" "$B" "$VD" "$DL" "$DT" "$coefL" "$coefT" "$(basename "$out")"
}

echo "Magboltz-derived pars at E = $EPOINT V/cm = $EFIELD_VM V/m, sampling $SAMPRATE MHz (scan: $SCAN/$PREFIX)"
emit 2.0  B20
emit 2.85 B285
emit 3.8  B38
echo
echo "NOTE: ATTPC.46Ar_3Hed_sim_B38.par is OVERWRITTEN -- the old one held the H2 placeholders."
