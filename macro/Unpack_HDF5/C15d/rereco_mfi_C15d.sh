#!/usr/bin/env bash
# C15d RE-RECO with MultiFit + TriplClust, dE/dx from the INTEGRATED pulse.
# Replaces rereco_hdb_C15d.sh and rereco_hdb2.
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
# ★★★ THE ALGORITHM CHOICE, MEASURED ON THE METRIC THAT MATTERS
#   The first attempt at this re-reco used MultiFit + HDBSCAN with space charge off, chosen on
#   TRACK and HIT counts: HDBSCAN finds +81 % tracks and 3.4x the hits-on-track, MultiFit only
#   +2 %/+4 %.  Those are the wrong numbers.  The quantity the analysis consumes is
#   PID-VALID TRACKS WHOSE VERTEX LANDS ON THE BEAM AXIS, because AtGenfitter gates on
#   AtSpyralResult::valid and AtSpyralPID rejects any vertex more than 25 mm off axis.  A config
#   can double the first two while losing a third of the last.  Measured one change at a time:
#
#     run_0039, 3,000 ev        tracks  hits/trk  valid  on-axis  on-axis 30-80 deg
#       max + tc  + SC           2,697     130    1,583   1,473        755   <- old production
#       max + tc  no SC          2,719     129      973     944          -
#       MF  + tc  + SC           2,694     138    1,822   1,723        825   <- BEST
#       MF  + hdb + SC           5,125     251    2,390   1,061        480
#       MF  + hdb no SC          5,128     251    2,366   1,064        483
#
#     run_0040, 6,000 ev        tracks  hits/trk  valid  on-axis
#       max + tc  + SC           5,391     126    3,122   2,925
#       MF  + tc  + SC           5,287     138    3,488   3,306   <- +13 % on-axis over production
#       MF  + tc  no SC          5,280     138    3,473   3,294   <- SC irrelevant WITH MultiFit
#       MF  + hdb + SC          10,324     240    4,706   2,179   <- -26 %
#
#   So: MULTIFIT IS THE REAL GAIN (+12-17 % usable tracks, +10 % hits/track) and HDBSCAN is a
#   REGRESSION here despite its track count.  Checked that HDBSCAN's deficit is not an artefact of
#   the selection: the Brho residual against the two-body elastic curve has the same shape in every
#   config and HDBSCAN is down in every bin, so it is depleted, not displaced.
#
#   SPACE CHARGE: the full 2x2, on-axis valid tracks, consistent across both runs --
#
#                      run_0039 (3,000 ev)   run_0040 (6,000 ev)
#       max + tc  SC on       1,473                2,925
#       max + tc  SC off        944                1,835     <- worth +56 % / +59 %
#       MF  + tc  SC on       1,723                3,306
#       MF  + tc  SC off      1,727                3,294     <- EXACTLY NEUTRAL
#
#   So SC was compensating for AtPSAMax's hit mis-placement; MultiFit resolves the overlapping
#   pulses properly and makes it redundant.  Default OFF, which also drops a foreign field map
#   (its LUTs are a1954's, and there is no a2091 one), 2.5x the output size and the CPU.
#   `$3=true` re-enables it.
#   ⚠ Do NOT conclude from this that SC is inert in general -- on a max-PSA reco it is not, and
#   track/hit counts cannot see it either way.
#
# ★★★ dE/dx FROM THE INTEGRATED PULSE (pidIntegral, $4, default true here).
#   ATTPCROOT measured dE/dx from the pulse PEAK; Spyral integrates. That single difference is
#   why the Spyral PID plane separates better. AtPSAMultiFit DID have a trace integral but built
#   it from the NOMINAL peaking time while the fit had measured another, making it exactly
#   2.9618x the amplitude for every one of 57,151 hits (rms/mean 0.0000) -- fixed 2026-11-05.
#   With the fix and a floated peaking time the integral is real and reproduces Spyral:
#       region     peak -> integral      Spyral vs ATTPCROOT
#       t               x1.27                   x1.20
#       alpha           x1.25                   x1.24
#   d-t resolving power 2.49 -> 4.31, deuteron WIDTH unchanged (0.1077 -> 0.1105): this buys
#   SEPARATION, not resolution. Costs +31 % CPU. The reconstruction itself barely moves
#   (2743 -> 2732 tracks, 153.6 -> 156.1 hits/track), so it is a charge definition, not a new reco.
#   ⚠ The integral is NON-LINEAR in true dE/dx at high ionisation: on the parameter-free 4He
#   spacing (ln(4/2)=0.693 at every Brho) the peak gives +0.66..0.70 and the integral +0.83..0.90.
#   Good for GATING; keep the peak-based plane in C15d_reco_mf for dE/dx used as a number.
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
OUT=/home/yassid/Data/C15d_reco_mfi
LOG=$OUT/logs
NPAR="${1:-16}"
DOSC="${3:-false}"        # space charge; see the table above -- irrelevant once MultiFit is used
PIDINT="${4:-true}"       # dE/dx from the integrated pulse; implies a floated peaking time
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
  root -b -q "$HERE/unpackReco_C15d.C(\"$r\",0,false,\"$OUT/\",\"$STAGE/\",\"ATTPC.C15d_a2091_D2.par\",\"ATTPC_D300torr_v2_geomanager.root\",20.0,15.0,7.5,true,false,\"\",30,0.0,2.85,0,\"multifit\",\"tc\",$DOSC,false,$PIDINT)" >"$LOG/${r}.log" 2>&1 || true
  local n; n=$(grep -aoE 'Hit last event at: [0-9]+' "$LOG/${r}.log" | tail -1 | grep -oE '[0-9]+$')
  if [ -s "$OUT/${r}_reco.root" ] && [ -n "$n" ]; then
    echo "[$(date +%H:%M:%S)] $r OK  $(du -h "$OUT/${r}_reco.root"|cut -f1)  ${n} ev  $(( (SECONDS-t0)/60 )) min"
  else
    echo "[$(date +%H:%M:%S)] $r FAILED  (no end-of-run marker; see $LOG/${r}.log)"
  fi
  rm -f "$STAGE/$r.h5"
}
# WARNING: DOSC MUST BE EXPORTED. Workers run `one` through `bash -c`, which does not inherit
# `set -u`, so an unexported DOSC expands to empty, cling sees a trailing comma, and every run
# dies in 5 s with "expected expression" while the driver correctly reports FAILED.
export -f one; export OUT STAGE LOG HERE DOSC PIDINT

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

echo "=== C15d re-reco MF+TC: ${#RUNS[@]} runs, $NPAR workers, cap $((CAP/1024/1024/1024)) GB, doSC=$DOSC, pidIntegral=$PIDINT -> $OUT ==="
echo "=== $(printf '%s ' "$(date)") already finished: $(for i in "${RUNS[@]}"; do done_already "$(printf run_%04d "$i")" && echo x; done | wc -l) ==="
stager & SPID=$!
trap 'kill $SPID 2>/dev/null' EXIT
printf '%s\n' "${RUNS[@]}" | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
wait $SPID 2>/dev/null || true
echo "=== DONE: $(for i in "${RUNS[@]}"; do done_already "$(printf run_%04d "$i")" && echo x; done | wc -l)/${#RUNS[@]} complete ==="
