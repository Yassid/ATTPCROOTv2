#!/usr/bin/env bash
# TC variant (2026-11-14): praType "tc" (radius 15, distance 7.5 = adopted unpackReco_a1975_UKF.C), else = reco_tau_h2_local.sh.
# reco_tau_h2.sh with LOCAL STAGING (2026-11-11). Same reco/PID arguments, bit-for-bit.
# Why: /mnt/f is an external HDD over 9p; 6 jobs reading 8 GB h5 files at once seek-thrash it
# (24 MB read in 8 min per job). Here ONE stager copies h5 files sequentially to the local ext4
# disk, NPAR workers reconstruct from the local copy into a local dir, and the outputs are moved
# back to /mnt/f/a1975/reco_h2_tau/ with the same .taudone/.done markers. Every /mnt/f bulk
# transfer (stager copy, output move) holds one flock, so the HDD only ever streams one file.
#
#   NPAR=8 ./reco_tau_h2_local.sh [runs...]    (runs that already have pid_width_XXXX.root.done are skipped)
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
set -uo pipefail

FH5=/mnt/f/a1975/h5/
REC=/mnt/f/a1975/reco_h2_tctau/
LOG=/mnt/f/a1975/logs_h2_tctau/
STG=/home/yassid/h2stage/
SH5=${STG}h5/
SREC=${STG}rec/
LOCK=${STG}f.lock
PAR=ATTPC.a1954.par
NPAR=${NPAR:-8}
MAXSTAGED=${MAXSTAGED:-2}   # 2026-11-12: 10 x ~8-21 GB filled C: (vhdx)
mkdir -p "$REC" "$LOG" "$SH5" "$SREC"

PEND=""
for r in ${*:-$(seq -f %04g 106 189)}; do
  [ -e "${REC}pid_width_${r}.root.done" ] && { echo "[have] run_$r"; continue; }
  PEND="$PEND $r"
done
echo "pending:$PEND"

stager() {
  for r in $PEND; do
    f="run_$r.h5"
    [ -e "${SH5}$f.ready" ] && continue
    while [ "$(ls "$SH5" | grep -c '\.ready$')" -ge "$MAXSTAGED" ]; do sleep 20; done
    # 2026-11-14: right after a reboot F: (drvfs) returned "Permission denied" on every h5 for a
    # while; a failed stage leaves no .ready and its worker waits forever -> retry for ~30 min.
    ok=0
    for try in 1 2 3 4 5 6; do
      flock "$LOCK" cp "${FH5}$f" "${SH5}$f.part" && mv "${SH5}$f.part" "${SH5}$f" && touch "${SH5}$f.ready" \
        && { ok=1; break; }
      rm -f "${SH5}$f.part"; echo "[retry stage] run_$r try $try"; sleep 300
    done
    [ $ok = 1 ] || echo "[FAIL stage] run_$r"
    echo "[staged] run_$r $(date '+%H:%M:%S')"
  done
}

one() {
  r="run_$1"
  until [ -e "${SH5}$r.h5.ready" ]; do sleep 10; done
  rc="${SREC}${r}_multifit_reco.root"
  pw="${SREC}pid_width_$1.root"
  echo "[reco ] $r $(date '+%H:%M:%S')"
  root -l -b -q "unpackReco_multifit.C(\"$r\",0,false,\"$SREC\",\"$SH5\",true,false,\"multifit\",0,20,\"tc\",20,8,-1,0.1,\"$PAR\",kTRUE)" \
      > "${LOG}reco_${r}.log" 2>&1
  # 2026-11-12: a reco KILLED mid-run leaves a big partial file and no segfault; only a normal exit
  # prints ROOT's "Real time" line (0172/0173 were marked done truncated that way).
  if grep -q 'segmentation violation' "${LOG}reco_${r}.log" || ! tail -5 "${LOG}reco_${r}.log" | grep -q 'Real time' \
     || [ ! -s "$rc" ] \
     || [ "$(stat -c%s "$rc" 2>/dev/null || echo 0)" -lt 20000000 ]; then
    # 2026-11-12: drop the local h5 COPY (original on /mnt/f) so a failed run cannot hold a
    # staging slot forever and fill C: (the ext4.vhdx); the partial reco + log stay for inspection.
    rm -f "${SH5}$r.h5" "${SH5}$r.h5.ready"
    echo "[FAIL reco] $r"; return 1
  fi
  rm -f "${SH5}$r.h5" "${SH5}$r.h5.ready"
  echo "[pid  ] $r $(date '+%H:%M:%S')"
  root -l -b -q "pid_width_dp.C(\"$rc\",\"$pw\")" > "${LOG}pid_${r}.log" 2>&1
  if [ ! -s "$pw" ]; then echo "[FAIL pid] $r"; return 1; fi
  flock "$LOCK" bash -c "cp '$rc' '${REC}.${r}.tmp' && mv '${REC}.${r}.tmp' '${REC}${r}_multifit_reco.root' && cp '$pw' '${REC}pid_width_$1.root'" \
    || { echo "[FAIL move] $r"; return 1; }
  touch "${REC}${r}_multifit_reco.root.taudone" "${REC}pid_width_$1.root.done"
  rm -f "$rc" "$pw"
  echo "[done ] $r $(date '+%H:%M:%S')  $(du -h "${REC}${r}_multifit_reco.root" | cut -f1)"
}
export -f one; export REC LOG PAR SH5 SREC LOCK

stager &
printf '%s\n' $PEND | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
wait
echo "=== reco_tau_h2 finished $(date '+%H:%M:%S') ==="
