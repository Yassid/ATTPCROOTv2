#!/usr/bin/env bash
# C15d RE-RECO with MultiFit + HDBSCAN, no space charge.  Replaces rereco_hdb_C15d.sh.
#
# WHY A REWRITE (measured 2026-11-03 on the first 12 runs of the v1 driver):
#   v1 processed a CHUNK of 6, waited for all 6, then started the next chunk.  Runs are wildly
#   unequal -- run_0013 has 2,070 events and run_0018 has 65,495 -- so the barrier made every
#   chunk as slow as its longest run.  With 6 workers on a 32-core box the CPU sat 84 % IDLE and
#   iowait was 3 %: not disk-bound, just under-parallelised.  Chunk 1 took 88 min for 6 runs,
#   i.e. ~25 h for 104 runs.  A queue with 16 workers projects ~5 h.
#
#   v1's success check was also wrong: it took the LAST "Event Number" line, which is
#   ATTPCROOT's 2^64-1 end-of-run sentinel, so every run printed
#     FAILED (events=18446744073709551615)
#   plus a bash "integer expression expected" error, while writing perfectly good output.
#   The real end-of-run marker is "[WARN] Hit last event at: N" -- used here, and it doubles as
#   the resume test, so a half-written .root from a killed worker is NOT mistaken for a finished
#   one (a size check alone cannot tell them apart).
#
# ★ STAGES TO NVMe, continuously rather than per chunk.  Raw runs are 13-30 GB and the whole
#   dataset is ~2 TB, which does not fit, so the stager runs ahead of the workers and blocks
#   while staged bytes exceed CAP.  Workers delete their own .h5 when done.  A worker that holds
#   a file is computing, never waiting, so the queue cannot deadlock on CAP.
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
RAW="/media/yassid/Seagate Hub/ATTPC/Data/a2091"
STAGE=/home/yassid/Data/C15d_stage
OUT=/home/yassid/Data/C15d_reco_hdb
LOG=$OUT/logs
NPAR="${1:-16}"
CAP=$(( ${2:-420} * 1024*1024*1024 ))        # max bytes resident in $STAGE
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
mkdir -p "$STAGE" "$OUT" "$LOG"; cd "$HERE"
mapfile -t RUNS < <(grep -vE '^\s*(#|$)' runs_a2091_have.txt | tr -s ' \n' '\n' | grep -oE '[0-9]+' | sed 's/^0*//')

# A run is finished iff its log carries the end-of-run marker AND the .root is non-trivial.
done_already(){ local r="$1"
  [ -s "$OUT/${r}_reco.root" ] && grep -aq 'Hit last event at' "$LOG/${r}.log" 2>/dev/null; }
export -f done_already

one(){ local i="$1"; local r; r=$(printf run_%04d "$i")
  done_already "$r" && { echo "[$(date +%H:%M:%S)] $r done already"; return 0; }
  local w=0
  while [ ! -s "$STAGE/$r.h5" ]; do
    [ -f "$STAGE/$r.missing" ] && { echo "[$(date +%H:%M:%S)] $r NO RAW FILE"; return 0; }
    sleep 10; w=$((w+10))
    [ $w -gt 10800 ] && { echo "[$(date +%H:%M:%S)] $r STAGE TIMEOUT"; return 0; }
  done
  rm -f "$OUT/${r}_reco.root"                 # drop any partial from an earlier kill
  local t0=$SECONDS
  root -b -q "$HERE/unpackReco_C15d.C(\"$r\",0,false,\"$OUT/\",\"$STAGE/\",\"ATTPC.C15d_a2091_D2.par\",\"ATTPC_D300torr_v2_geomanager.root\",20.0,15.0,7.5,true,false,\"\",30,0.0,2.85,0,\"multifit\",\"hdbscan\",false)" >"$LOG/${r}.log" 2>&1 || true
  local n; n=$(grep -aoE 'Hit last event at: [0-9]+' "$LOG/${r}.log" | tail -1 | grep -oE '[0-9]+$')
  if [ -s "$OUT/${r}_reco.root" ] && [ -n "$n" ]; then
    echo "[$(date +%H:%M:%S)] $r OK  $(du -h "$OUT/${r}_reco.root"|cut -f1)  ${n} ev  $(( (SECONDS-t0)/60 )) min"
  else
    echo "[$(date +%H:%M:%S)] $r FAILED  (no end-of-run marker; see $LOG/${r}.log)"
  fi
  rm -f "$STAGE/$r.h5"
}
export -f one; export OUT STAGE LOG HERE

stager(){ local i r sz used
  for i in "${RUNS[@]}"; do
    r=$(printf run_%04d "$i")
    done_already "$r" && continue
    [ -s "$STAGE/$r.h5" ] && continue
    [ -f "$RAW/$r.h5" ] || { touch "$STAGE/$r.missing"; continue; }
    sz=$(stat -c%s "$RAW/$r.h5")
    while :; do
      used=$(du -sb "$STAGE" 2>/dev/null | cut -f1)
      [ $((used + sz)) -lt $CAP ] && break
      sleep 20
    done
    cp "$RAW/$r.h5" "$STAGE/$r.h5.part" && mv "$STAGE/$r.h5.part" "$STAGE/$r.h5"
  done
  echo "[$(date +%H:%M:%S)] stager: all runs staged or skipped"
}

echo "=== C15d re-reco v2: ${#RUNS[@]} runs, $NPAR workers, stage cap $((CAP/1024/1024/1024)) GB -> $OUT ==="
echo "=== $(printf '%s ' "$(date)") already finished: $(for i in "${RUNS[@]}"; do done_already "$(printf run_%04d "$i")" && echo x; done | wc -l) ==="
stager & SPID=$!
trap 'kill $SPID 2>/dev/null' EXIT
printf '%s\n' "${RUNS[@]}" | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
wait $SPID 2>/dev/null || true
echo "=== DONE: $(for i in "${RUNS[@]}"; do done_already "$(printf run_%04d "$i")" && echo x; done | wc -l)/${#RUNS[@]} complete ==="
