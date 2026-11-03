#!/usr/bin/env bash
# ★ SUPERSEDED by rereco_hdb2_C15d.sh -- this version serialises on a chunk barrier (6 workers
#   on a 32-core box, CPU 84 % idle, ~25 h for 104 runs) and its success check reads the
#   2^64-1 end-of-run sentinel, so it prints FAILED for every run while writing good output.
#   Kept only for provenance: the 2026-11-03 runs 0013-0024 were produced with it.
# C15d RE-RECO with HDBSCAN. The 2026-09-03 production used AtPSAMax + AtTrackFinderTC + space
# charge; measured on 300 events of run_0017 that configuration finds 278 tracks carrying 38,218
# hits, against 530 / 131,557 with MultiFit + HDBSCAN.
#
#   ISOLATED (same 300 events, one change at a time):
#     max+tc+SC    278 tracks   38,218 hits-on-track   122/track   <- production
#     MF +tc+SC    283 (+2%)    39,898 (+4%)           127         <- MultiFit alone: marginal
#     max+hdb+SC   504 (+81%)  122,802 (+221%)         157         <- HDBSCAN alone: the whole gain
#     MF +hdb+SC   525         131,970                 153
#     MF +hdb noSC 530         131,557                 152         <- SC changes nothing
#
# So HDBSCAN is the reason, not the pulse analysis. Space charge is dropped because it measurably
# does nothing here AND its LUTs are the a1954 experiment's field map, not this one's.
#
# ★ STAGES TO NVMe. Spyral on this same drive measured 656 kB/s with 8 workers because its reads
# are small and random; ATTPCROOT streams better but still ran 122 s wall for 4.6 s CPU on a cold
# Seagate file. Sequential copy is 200 MB/s, so runs are copied first and workers read local disk.
set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
RAW="/media/yassid/Seagate Hub/ATTPC/Data/a2091"
STAGE=/home/yassid/Data/C15d_stage
OUT=/home/yassid/Data/C15d_reco_hdb
LOG=$OUT/logs
CHUNK="${1:-6}"; NPAR="${2:-6}"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
mkdir -p "$STAGE" "$OUT" "$LOG"; cd "$HERE"
mapfile -t RUNS < <(grep -vE '^\s*(#|$)' runs_a2091_have.txt | tr -s ' \n' '\n' | grep -oE '[0-9]+' | sed 's/^0*//')
echo "=== C15d re-reco: ${#RUNS[@]} runs, chunks of $CHUNK, $NPAR workers -> $OUT ==="
one(){ local i="$1"; local r=$(printf run_%04d "$i")
  [ -s "$OUT/${r}_reco.root" ] && [ $(stat -c%s "$OUT/${r}_reco.root") -gt 1000000 ] && { echo "[$(date +%H:%M:%S)] $r exists"; return 0; }
  [ -f "$STAGE/$r.h5" ] || { echo "[$(date +%H:%M:%S)] $r NOT STAGED"; return 0; }
  local t0=$SECONDS
  root -b -q "$HERE/unpackReco_C15d.C(\"$r\",0,false,\"$OUT/\",\"$STAGE/\",\"ATTPC.C15d_a2091_D2.par\",\"ATTPC_D300torr_v2_geomanager.root\",20.0,15.0,7.5,true,false,\"\",30,0.0,2.85,0,\"multifit\",\"hdbscan\",false)" >"$LOG/${r}.log" 2>&1
  local n=$(grep -aoE "Event Number : +[0-9]+" "$LOG/${r}.log" | tail -1 | grep -oE "[0-9]+$")
  if [ -s "$OUT/${r}_reco.root" ] && [ "${n:-0}" -ge 100 ]; then
     echo "[$(date +%H:%M:%S)] $r OK $(du -h "$OUT/${r}_reco.root"|cut -f1) ${n} ev in $(( (SECONDS-t0)/60 )) min"
  else echo "[$(date +%H:%M:%S)] $r FAILED (events=${n:-none})"; fi
}
export -f one; export OUT STAGE LOG HERE
stage_c(){ for i in "$@"; do f=$(printf run_%04d "$i").h5; [ -s "$STAGE/$f" ] && continue
             [ -f "$RAW/$f" ] || continue; cp "$RAW/$f" "$STAGE/$f.part" && mv "$STAGE/$f.part" "$STAGE/$f"; done; }
idx=0; pf=""; cn=0
while [ $idx -lt ${#RUNS[@]} ]; do
  cn=$((cn+1)); cur=("${RUNS[@]:$idx:$CHUNK}"); nxt=("${RUNS[@]:$((idx+CHUNK)):$CHUNK}")
  [ -n "$pf" ] && wait "$pf" 2>/dev/null || true
  echo "[$(date +%H:%M:%S)] chunk $cn: staging ${cur[*]}"; stage_c "${cur[@]}"
  [ ${#nxt[@]} -gt 0 ] && { stage_c "${nxt[@]}" & pf=$!; } || pf=""
  printf '%s\n' "${cur[@]}" | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
  for i in "${cur[@]}"; do rm -f "$STAGE/$(printf run_%04d $i).h5"; done
  echo "[$(date +%H:%M:%S)] chunk $cn done -- $((idx+${#cur[@]}))/${#RUNS[@]}; free $(df -h /home/yassid/Data|tail -1|awk '{print $4}')"
  idx=$((idx+CHUNK))
done
[ -n "$pf" ] && wait "$pf" 2>/dev/null || true
echo "=== DONE: $(ls -1 $OUT/*_reco.root 2>/dev/null|wc -l)/${#RUNS[@]} ==="
