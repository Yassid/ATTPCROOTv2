#!/usr/bin/env bash
# Fit the SAME events at several truncation percentages, each into its own directory.
# 0 is the NULL TEST: it must reproduce the production fit exactly, and if it does not, nothing
# else here means anything. Judge the rest on Ex, never on chi2 -- chi2 falls monotonically as
# clusters are removed, which is how "best chi2" once returned a 13%-biased KE.
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export GENFIT=/home/yassid/fair_install/GenFit
export LD_LIBRARY_PATH=$GENFIT/lib:${LD_LIBRARY_PATH:-}
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
set -uo pipefail
RUN="${RUN:-run_0020_multifit}"; NEV="${NEV:-3000}"
one() {
  pct="$1"; out="/mnt/f/a1975/gf_dp_trunc/p${pct}/"; mkdir -p "$out"
  root -l -b -q "fitGenfitter_a1975_deuterium.C(\"$RUN\",$NEV,\"/mnt/f/a1975/reco_d2_dv1104/\",\"\",\"$out\",\
-2.85,2,5,\"pid/proton_d2_dv1104.json\",4.0,10.0,170.0,kTRUE,kTRUE,2212,1.00782503207,1,\"p\",\"_reco\",\
\"ATTPC_D300torr_v2_geomanager.root\",kTRUE,0,2,\"ATTPC.a1975_deuterium_dv1104.par\",kFALSE,kFALSE,kFALSE,\
\"proton_D2_300torr.txt\",kTRUE,kTRUE,kTRUE,kFALSE,${pct})" > "/mnt/f/a1975/gf_dp_trunc/fit_p${pct}.log" 2>&1
  echo "[done] pct=$pct  $(date '+%H:%M:%S')"
}
export -f one; export RUN NEV
printf '%s\n' ${*:-25 50 75 90} | xargs -P 4 -I{} bash -c 'one "$@"' _ {}
echo "=== trunc scan finished $(date '+%H:%M:%S') ==="
