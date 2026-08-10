#!/usr/bin/env bash
# Gain scan for the a1975 (p,p) simulation, against one number the data measures directly.
#
# WHY. The simulated acceptance collapsed at forward theta_cm, where the recoil proton carries
# 0.5-2 MeV -- but the DATA has 138,000 events there with the same lab kinematics, so the detector
# clearly reconstructs them. The discrepancy is in the simulation:
#
#     fraction of fitted tracks with KE < 2 MeV     data 57.8 %     simulation 7.0 %
#
# a factor of eight. Gain governs whether a short, low-energy track clears threshold at all, and
# the gain in the par file is inherited from a1954 and has never been checked against a1975. The
# same mistake there had the gain 15x too low.
#
# The scan reruns the SAME generated events through reco + fit at each gain, so nothing but the
# digitisation changes and the comparison is like-for-like.
#
# The target is the data's low-energy fraction. It is a necessary condition, not a sufficient one:
# a gain that reproduces it could still get the track width wrong, which is what decided the a1954
# choice (2.5e6 matched the charge but flattened the width's sqrt(z) growth, so 150k was kept).
# Track width per gain is printed alongside for that reason.
#
#   ./gain_scan_pp.sh [simFile] [outDir]
set -eo pipefail
SIM=${1:-/mnt/f/a1975_C16_pp_sim_el300.root}
OUT=${2:-/mnt/f/a1975_C16_pp_gain}
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
HERE=$REPO/macro/Simulation/ATTPC/16C_pp_a1975
mkdir -p "$OUT"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"

for G in 150000 400000 1000000 2500000 5000000; do
  PAR="ATTPC.a1975_C16_sim.par"
  [ "$G" != "150000" ] && PAR="ATTPC.a1975_C16_sim_g$G.par"
  J="g$G"
  if [ -f "$OUT/$J.marker" ]; then echo "[$(date +%H:%M:%S)] $J already done"; continue; fi
  echo "[$(date +%H:%M:%S)] gain $G : reco"
  cd "$REPO/macro/Simulation/ATTPC/16C_pd"
  root -b -q -l "run_reco_C16pd.C(\"$SIM\",\"$OUT/${J}_reco.root\",\"$PAR\",20,20,8,0,20.0,\"mover\")" \
       > "$OUT/${J}_reco.log" 2>&1
  grep -q "sim reco done" "$OUT/${J}_reco.log" || { echo "  $J RECO_FAILED"; continue; }
  echo "[$(date +%H:%M:%S)] gain $G : fit"
  cd "$REPO/macro/Unpack_HDF5/a1975/UKF"
  root -b -q -l "pipeline/fitGenfitter_a1975.C(\"$J\",-1,\"$OUT/\",\"\",\"$OUT/\",-2.85,2,5,\"proton_H2_catima.txt\",\"\",4.0,kFALSE,2212,1.00727647,1,kFALSE,\"ATTPC_H300torr_RT\")" \
       > "$OUT/${J}_fit.log" 2>&1
  # a written file is not a finished job
  [ -s "$OUT/${J}_genfitter.root" ] && echo COMPLETED > "$OUT/$J.marker" || echo "  $J FIT_FAILED"
  rm -f "$OUT/${J}_reco.root"   # ~1.4 GB each and nothing downstream needs it once the fit exists
done
echo "[$(date +%H:%M:%S)] SCANDONE"
