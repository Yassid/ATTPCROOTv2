#!/bin/bash
# Scan UKF σ_p/p vs pion momentum for HYDRA Prototype.
# For each p in PLIST: run sim → digi → UKF → record σ_p/p.
#
# Usage: ./scan_p.sh   (run from macro/Simulation/ATTPC/HYDRA)

set -e

# Source the ATTPCROOT build config so ROOT, VMCWORKDIR, etc. are in path.
# Self-contained: works whether or not the parent shell already sourced it.
ATTPCROOT_BUILD="/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build"
if [ -f "$ATTPCROOT_BUILD/config.sh" ]; then
   # shellcheck disable=SC1091
   source "$ATTPCROOT_BUILD/config.sh"
fi

# Pion mass [MeV/c²]
MASS=139.57039
NEVT=2000

# Momentum points [MeV/c]
PLIST=(200 400 600 800 1000 1200)

OUT=data/scan_p_results.csv
echo "p_MeV,KE_MeV,n_fit,n_thr,mean_pMC,bias,sigma_p_over_p" > "$OUT"

for P in "${PLIST[@]}"; do
   # KE = sqrt(p² + m²) − m
   KE=$(python3 -c "import math; m=$MASS; p=$P; print(f'{math.sqrt(p*p + m*m) - m:.2f}')")
   KEMIN=$(python3 -c "print(f'{$KE - 5:.2f}')")
   KEMAX=$(python3 -c "print(f'{$KE + 5:.2f}')")
   SUF="_p${P}"

   echo "===== p=$P MeV/c, KE=$KE MeV (window $KEMIN-$KEMAX) ====="

   rm -f data/HYDRAsim${SUF}.root data/HYDRApar${SUF}.root data/output_digi${SUF}.root data/output_ukf_HYDRA${SUF}.root

   echo "[sim]"
   root -b -q "HYDRA_sim.C($NEVT, -1, $KEMIN, $KEMAX, \"TGeant4\", 2.0, \"$SUF\", 1)" > /tmp/scan_sim${SUF}.log 2>&1
   echo "[digi]"
   root -b -q "run_digi_HYDRA.C($NEVT, 8.0, false, \"$SUF\", true, 6.0, \"$SUF\")" > /tmp/scan_digi${SUF}.log 2>&1
   echo "[ukf]"
   root -b -q "run_ukf_HYDRA.C($NEVT, -1, 0.1, 3, \"$SUF\", 1.0, \"$SUF\", 2.0)" > /tmp/scan_ukf${SUF}.log 2>&1

   # Extract σ_p/p from p_resolution analysis macro
   LINE=$(root -b -q "analysis/p_resolution.C(\"data/output_ukf_HYDRA${SUF}.root\",\"data/HYDRAsim${SUF}.root\",\"p=${P}\")" 2>&1 \
          | awk '/80-100/ {print}')
   # Format: "  80-100  <Nthr> <Nfit>  <pMC>  <bias>  <sigma> ( ... %)"
   NTHR=$(echo "$LINE"  | awk '{print $2}')
   NFIT=$(echo "$LINE"  | awk '{print $3}')
   PMC=$(echo "$LINE"   | awk '{print $4}')
   BIAS=$(echo "$LINE"  | awk '{print $5}')
   SIG=$(echo "$LINE"   | awk '{print $6}')

   echo "p=$P MeV/c → σ_p/p=$SIG  (N=$NFIT/$NTHR)"
   echo "$P,$KE,$NFIT,$NTHR,$PMC,$BIAS,$SIG" >> "$OUT"
done

echo "===== Done ====="
column -t -s, "$OUT"
