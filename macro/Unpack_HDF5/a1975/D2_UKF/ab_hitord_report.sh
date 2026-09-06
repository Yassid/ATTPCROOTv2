#!/usr/bin/env bash
# Waits for reco_hitord_ab.sh to finish each run, then runs the three A/B measurements on it.
# Writes one report per run to /mnt/f/a1975/logs_hitord/ab_report_<run>.txt.
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
set -uo pipefail
REC=/mnt/f/a1975/reco_d2_dv1104_hitord/
LOG=/mnt/f/a1975/logs_hitord/

for n in 0020 0016; do
  r="run_$n"; f="${REC}${r}_multifit_reco.root"; out="${LOG}ab_report_${r}.txt"
  # wait for the driver's completion marker, but give up if the reco died (no root process left)
  while [ ! -e "${f}.hitorddone" ]; do
    pgrep -f "unpackReco_multifit.*${r}" >/dev/null || { echo "[abort] $r: reco not running and no marker" > "$out"; break; }
    sleep 60
  done
  [ -e "${f}.hitorddone" ] || continue
  {
    echo "=========== A/B hit-ordering report: $r   $(date '+%F %H:%M:%S')"
    echo
    echo "--- 1. did the fix take effect? re-ordering the NEW file's hits should now be a no-op ---"
    root -l -b -q "hit_path_ratio.C(\"$f\",400,30)" 2>&1 | tail -6
    root -l -b -q "theta_order_test.C(\"$f\",400,30)" 2>&1 | tail -6
    echo
    echo "--- 2. control, for reference: the SAME measurements on the production reco ---"
    root -l -b -q "hit_path_ratio.C(\"/mnt/f/a1975/reco_d2_dv1104/${r}_multifit_reco.root\",400,30)" 2>&1 | tail -6
    echo
    echo "--- 3. the A/B itself: production vs re-reco, matched on (entry, trackID) ---"
    root -l -b -q "ab_hitord_compare.C(\"$r\")" 2>&1 | tail -14
  } > "$out" 2>&1
  echo "[report] $out"
done
echo "=== ab_hitord_report finished $(date '+%F %H:%M:%S') ==="
