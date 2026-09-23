#!/usr/bin/env bash
# C15p: FULL reconstruction from raw HDF5 with AtPSAMultiFit + HDBSCAN, space charge OFF.
#
#   ./reco_batch_C15p.sh [nparallel] [runlist]
#
# CONFIG (differs from every earlier a2091 reco):
#   PSA   AtPSAMultiFit, threshold 20     -- deconvolves overlapping pulses on a pad
#   CLEAN AtDirDeDxCleaner               -- noise removal before pattern recognition
#   PRA   AtTrackFinderHDBSCAN, "mover"   -- minClusterSize 20, minSamples 8
#   SC    OFF                             -- the a1954 distortion map is a different experiment's
#   PID   AtPIDTask in-chain              -- persists AtPIDEvent, so there is NO separate PID pass
#   geom  ATTPC_H300torr_RT, par ATTPC.a2091_C15.par
#
# ############################################################################################
# ★ WHY THIS STAGES TO NVMe INSTEAD OF READING THE SEAGATE DIRECTLY.
# The first attempt ran 20 workers straight off the USB drive and collapsed: 19 of 20 processes
# sat in D state, sda pinned at 99.7 % util delivering 3.1 MB/s with r_await 90 ms, and the
# projected wall time was 75 hours.
#
# The cause is NOT bandwidth. A single worker needs about 1 MB/s (measured: 300 events of
# run_0138 = 162 MB read in 158 s -- the job is CPU-bound in the MultiFit pulse fitting), so 20
# workers want only ~25 MB/s and the drive streams 62 MB/s on one thread. What kills it is the
# ACCESS PATTERN: 20 interleaved streams across 20 files spread over a 7 TB platter turn every
# read into a seek. Lowering the worker count fixes the disk and wastes the CPU -- at 3 workers
# the same job takes 60+ hours.
#
# So the raw file is copied to NVMe first, one at a time (sequential, the one thing the drive is
# good at), and the workers read only from NVMe, where interleaving is free. A background stager
# keeps a bounded queue ahead of the workers so copying overlaps with reconstruction, and each
# worker deletes its staged copy when the run is done, so the staging area never grows.
# ############################################################################################
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
NPAR="${1:-24}"
RUNLIST="${2:-$HERE/runs_pp_hdb.txt}"
SRC="/media/yassid/Seagate Hub/ATTPC/Data/a2091"
STAGE="${C15P_STAGE:-/home/yassid/Data/a2091_stage}"
OUT="${C15P_RECO_HDB:-/home/yassid/Data/a2091_C15p_reco_hdb}"
LOG="$OUT/logs"
MACRO="$REPO/macro/Unpack_HDF5/a2091/UKF/pipeline/unpackReco_C15.C"
# Ceiling on the staging area. Average run is ~20 GB, so this holds ~15 runs ahead of the
# workers. NVMe has ~924 GB free and the reco output is ~87 GB, so 450 GB still leaves headroom,
# and it keeps ~22 runs queued so the workers are never starved by the copy.
STAGE_MAX_GB="${STAGE_MAX_GB:-450}"

set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
[[ -n "${VMCWORKDIR:-}" ]] || { echo "ERROR: config.sh did not set VMCWORKDIR" >&2; exit 1; }
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
[[ -d "$SRC" ]] || { echo "ERROR: raw HDF5 dir not mounted: $SRC" >&2; exit 1; }
[[ -f "$MACRO" ]] || { echo "ERROR: missing $MACRO" >&2; exit 1; }
mkdir -p "$OUT" "$LOG" "$STAGE"
cd "$HERE"

mapfile -t RUNS < <(grep -vE '^\s*(#|$)' "$RUNLIST" | tr -s ' \n' '\n' | grep '^run_')
echo "=== C15p reco: ${#RUNS[@]} runs, $NPAR workers, staging via $STAGE (cap ${STAGE_MAX_GB} GB) ==="
echo "    PSA multifit(thr20) + AtDirDeDxCleaner + HDBSCAN(mover) + AtPIDTask, doSC=false"

done_already(){ [[ -s "$OUT/${1}_reco.root" ]] && [[ $(stat -c%s "$OUT/${1}_reco.root") -gt 1000000 ]]; }

# ---- stager: sequential copies, bounded queue, in worker order -----------------------------
stager(){
  for r in "${RUNS[@]}"; do
    done_already "$r" && continue
    [[ -f "$SRC/$r.h5" ]] || { echo "[stage] $r NO RAW"; touch "$STAGE/$r.norawX"; continue; }
    [[ -f "$STAGE/$r.ready" ]] && continue
    # wait for room
    while :; do
      used=$(du -sBG --apparent-size "$STAGE" 2>/dev/null | cut -dG -f1); used=${used:-0}
      [[ "$used" -lt "$STAGE_MAX_GB" ]] && break
      sleep 20
    done
    echo "[stage $(date +%H:%M:%S)] copying $r ($(du -h "$SRC/$r.h5" | cut -f1))"
    if cp "$SRC/$r.h5" "$STAGE/$r.h5.part" && mv "$STAGE/$r.h5.part" "$STAGE/$r.h5"; then
      touch "$STAGE/$r.ready"
    else
      echo "[stage] $r COPY FAILED"; rm -f "$STAGE/$r.h5.part"; touch "$STAGE/$r.norawX"
    fi
  done
  touch "$STAGE/.staging_complete"
  echo "[stage $(date +%H:%M:%S)] all copies done"
}

one(){
  local r="$1"
  done_already "$r" && { echo "[$(date +%H:%M:%S)] $r  exists, skipping"; return 0; }
  # wait for the stager (or for it to declare the run unavailable)
  local waited=0
  while [[ ! -f "$STAGE/$r.ready" ]]; do
    [[ -f "$STAGE/$r.norawX" ]] && { echo "[$(date +%H:%M:%S)] $r  unavailable, skipping"; return 0; }
    sleep 10; waited=$((waited+10))
    [[ $waited -gt 21600 ]] && { echo "[$(date +%H:%M:%S)] $r  TIMED OUT waiting for stage"; return 0; }
  done
  local t0=$SECONDS
  root -b -q -l "$MACRO(\"$r\",0,false,\"$OUT/\",\"$STAGE/\",false,false,\"multifit\",0,20,\"hdbscan\")" \
        >"$LOG/${r}_reco.log" 2>&1
  # A non-empty output file is NOT proof of success. run_0180 unpacks 1 event and writes a valid
  # 47 kB ROOT file when its HDF5 event count comes back as zero, and the old analysis carried
  # exactly that as a finished run (11 kB reco, 7.7 kB kin) without anyone noticing. So check the
  # EVENT COUNT the macro reported, not the file.
  local nev
  nev=$(grep -oE "Reconstructing [0-9]+ events" "$LOG/${r}_reco.log" 2>/dev/null | grep -oE "[0-9]+" | head -1)
  if [[ -s "$OUT/${r}_reco.root" ]] && [[ "${nev:-0}" -ge 100 ]]; then
     echo "[$(date +%H:%M:%S)] $r  OK  $(du -h "$OUT/${r}_reco.root" | cut -f1)  ${nev} events  in $(( (SECONDS-t0)/60 )) min"
  else
     echo "[$(date +%H:%M:%S)] $r  FAILED  (events=${nev:-none}, see $LOG/${r}_reco.log)"
  fi
  rm -f "$STAGE/$r.h5" "$STAGE/$r.ready"   # free the staging slot immediately
}
export -f one done_already; export OUT STAGE LOG MACRO

stager > "$LOG/stager.log" 2>&1 &
STAGER_PID=$!
trap 'kill $STAGER_PID 2>/dev/null || true' EXIT

printf '%s\n' "${RUNS[@]}" | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
wait $STAGER_PID 2>/dev/null || true
rm -f "$STAGE"/*.norawX "$STAGE/.staging_complete"
echo "=== reco done: $(ls -1 "$OUT"/*_reco.root 2>/dev/null | wc -l) / ${#RUNS[@]} files ==="
