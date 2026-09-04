#!/usr/bin/env bash
# Chunked reconstruction: ONE run at a time, split across many cores.
#
#   ./reco_chunked_C15d.sh [chunks] [runlist]
#
# WHY THIS EXISTS. The obvious way to use 32 cores is to reconstruct 32 runs at once. On this
# machine that is the WORST thing to do: the raw data is on a single spinning USB disk and each
# run is a 17-27 GB read, so concurrency turns sequential reads into seeks --
#
#     2 workers  -> 23.6 MB/s, 118 kB per I/O        12 workers -> 4.0 MB/s, 16 kB per I/O
#
# and per-worker CPU falls to 16 %, i.e. every worker waiting on the disk.
#
# But the drive streams at 253 MB/s sequentially, and a run whose file is in page cache
# reconstructs at 99 % CPU (measured, against 40 % cold). So the fast route is the opposite of the
# obvious one: read ONE run sequentially into cache, then let many cores share that warm file.
#
#   1  prefetch the run into page cache with a single sequential read
#   2  reconstruct it as N chunks in parallel, all served from RAM
#   3  hadd the chunks into <run>_reco.root
#
# Chunking is exact: chunk 0 + chunk 1 reproduce a single pass byte for byte (1200 events, 1085
# tracks, 132451 hits, identical sum of hit z), and hadd preserves that.
#
# ⚠ The chunk offset goes to AtUnpacker::SetInitialEventID BEFORE the unpacker is moved into the
# task. FairRunAna::Run(first, last) does NOT drive the HDF5 unpacker -- passing a start event
# there is silently ignored and every chunk re-reads the same events.
#
# ⚠ MEMORY. Each chunk is ~1.1 GB resident and the warm file needs 17-27 GB of page cache, against
# 62 GB total. 12 chunks is about the limit; more evicts the very file we prefetched.

set -eo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"
NCHUNK="${1:-12}"
RUNLIST="${2:-$HERE/runs_a2091_d2.txt}"
RAW="/media/yassid/Seagate Hub/ATTPC/Data/a2091/"
RECO="${C15D_RECO:-/home/yassid/C15d_reco}"
LOG="${C15D_LOGS:-/home/yassid/C15d_logs}"
PAR="ATTPC.C15d_a2091_D2.par"
GEO="ATTPC_D300torr_v2_geomanager.root"
MIN_FREE_GB=40
H5LS="/home/yassid/fair_install/hdf5-1.10.4-inst/bin/h5ls"
MAX_SANE_EVENTS=2000000   # a real AT-TPC run is 1e3-1e5 events; anything above this is corrupt metadata

# ★ HDF5 FILE LOCKING OFF. The raw data lives on an exFAT volume, which has no real file
# locking, and 12 chunks open the SAME .h5 read-only at the same moment. HDF5 1.10.4 then
# intermittently fails with "bad symbol table node signature" on /meta and the unpacker derives a
# garbage event range and calls LOG(fatal). It is a RACE, not corruption: h5ls reads /meta of the
# very same files without complaint. It cost 12 of 104 runs on the first full pass.
export HDF5_USE_FILE_LOCKING=FALSE
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
[[ -n "${VMCWORKDIR:-}" ]] || { echo "ERROR: config.sh did not set VMCWORKDIR" >&2; exit 1; }
mkdir -p "$RECO" "$LOG" "$RECO/.chunk"
cd "$HERE"

mapfile -t RUNS < <(grep -vE '^\s*(#|$)' "$RUNLIST" | tr -s ' \n' '\n' | grep '^run_')
echo "=== chunked reco: ${#RUNS[@]} runs, $NCHUNK chunks per run ==="
echo "  raw  : $RAW"
echo "  out  : $RECO"
echo "  par  : $PAR"

for run in "${RUNS[@]}"; do
   out="$RECO/${run}_reco.root"
   [[ -s "$out" ]] && { echo "[$run] exists, skipping"; continue; }
   raw="$RAW$run.h5"
   [[ -s "$raw" ]] || { echo "[$run] no raw file, skipping"; continue; }

   free_gb=$(df -BG --output=avail "$RECO" | tail -1 | tr -d ' G')
   (( free_gb < MIN_FREE_GB )) && { echo "[$run] SKIP: only ${free_gb} GB free"; continue; }

   t0=$SECONDS
   # --- 1. how many events? probe once; the unpacker prints it during Init -------------------
   # ★ THE `|| true` IS LOAD-BEARING. With `set -eo pipefail`, a probe that prints no
   # "Run contains" line makes grep return 1, so the pipeline returns 1, so the ASSIGNMENT
   # returns 1, and set -e kills the whole script -- before the `[[ -n "$nev" ]]` check below
   # that exists to handle exactly that case. The guard was unreachable. That is how a 104-run
   # pass silently stopped after 39 runs with no error message: run_0059's probe came back empty
   # and the script simply exited, and the caller's `|| true` made it look like normal completion.
   nev=$(root -b -q "$HERE/unpackReco_C15d.C(\"$run\",1,false,\"$RECO/.chunk/probe_\",\"$RAW\",\"$PAR\",\"$GEO\")" 2>/dev/null \
         | grep -oE "Run contains [0-9]+" | grep -oE "[0-9]+" | head -1) || true
   rm -f "$RECO/.chunk/probe_${run}_reco.root"
   if [[ -z "$nev" || "$nev" -le 0 ]]; then
      # The probe printed no event count at all. Try the same /get recovery used for corrupt
      # metadata before giving up -- two datasets per event (evt<N>_data, evt<N>_header).
      echo "[$run] probe gave no event count -- counting /get datasets instead"
      real=$(timeout 600 "$H5LS" "$raw/get" 2>/dev/null | wc -l) || true
      real=$(( ${real:-0} / 2 ))
      if (( real > 0 && real <= MAX_SANE_EVENTS )); then
         echo "[$run] recovered: $real events from the dataset count"
         nev=$real
      else
         echo "[$run] PROBE FAILED and unrecoverable -- skipping, continuing with the rest"
         echo "$run" >> "$HERE/runs_probe_failed.txt"
         continue
      fi
   fi
   # ★ SANITY-CHECK THE EVENT COUNT. Some files have corrupt metadata and the unpacker reports it
   # verbatim: run_0019 claims 127,979,076,285,538 events. Chunking on that asks for datasets like
   # evt95984350326776_header, which do not exist, and HDF5 prints a full multi-line error stack
   # PER ATTEMPT in a loop that never terminates -- 127 GB per chunk log, 1.4 TB over 12 chunks,
   # which filled the disk and hung every worker in D state. A real AT-TPC run is at most a few
   # hundred thousand events.
   if (( nev > MAX_SANE_EVENTS )); then
      # RECOVERABLE: the DATA is fine, only /meta is wrong. /get holds two datasets per event
      # (evtN_data and evtN_header), so counting its children and halving gives the true count.
      # run_0019: metadata claims 1.3e14, /get has 14568 children = 7284 real events.
      echo "[$run] bad metadata ($nev events) -- counting /get datasets instead"
      real=$(timeout 600 "$H5LS" "$raw/get" 2>/dev/null | wc -l) || true
      real=$(( real / 2 ))
      if (( real > 0 && real <= MAX_SANE_EVENTS )); then
         echo "[$run] recovered: $real events from the dataset count"
         nev=$real
      else
         echo "[$run] UNRECOVERABLE: /get gave $real events -- skipping"
         echo "$run" >> "$HERE/runs_bad_metadata.txt"
         continue
      fi
   fi

   # --- 2. prefetch the whole file into page cache, sequentially ------------------------------
   dd if="$raw" of=/dev/null bs=8M 2>/dev/null || true
   t_pf=$((SECONDS - t0))

   # --- 3. N chunks in parallel, all served from RAM ------------------------------------------
   per=$(( (nev + NCHUNK - 1) / NCHUNK ))
   pids=()
   for ((c=0; c<NCHUNK; c++)); do
      first=$(( c * per ))
      (( first >= nev )) && break
      # ★ CLAMP THE LAST CHUNK. per = ceil(nev/NCHUNK), so NCHUNK*per OVERRUNS nev -- for
      # run_0110, 12 x 2078 = 24936 against 24929 real events. The merged file then carries a few
      # phantom events at the end, its length no longer matches the IC file, and
      # make_points_C15d.C refuses the join: "IC has 24930 entries vs 24936 reco events". That
      # silently costs the run its beam gate. 22 of 35 runs were affected.
      this=$per
      (( first + this > nev )) && this=$(( nev - first ))
      # head -c caps each chunk log at 50 MB. Belt and braces against the runaway above: even
      # with the event-count guard, any per-event error loop would otherwise fill the disk.
      root -b -q "$HERE/unpackReco_C15d.C(\"$run\",$this,false,\"$RECO/.chunk/$(printf c%02d_ "$c")\",\"$RAW\",\"$PAR\",\"$GEO\",20.0,15.0,7.5,kTRUE,kFALSE,\"\",30,0.0,2.85,$first)" 2>&1 \
         | head -c 50000000 >"$LOG/${run}_$(printf c%02d "$c").log" &
      pids+=($!)
      sleep 0.4   # stagger the opens; simultaneous first-access is what trips the metadata cache
   done
   fail=0
   for p in "${pids[@]}"; do wait "$p" || fail=1; done

   # --- 4. merge ------------------------------------------------------------------------------
   # ★ MERGE IN NUMERIC ORDER, EXPLICITLY. A glob of c*_ expands ALPHABETICALLY --
   # c0_, c10_, c11_, c1_, c2_ ... -- so hadd concatenated the chunks out of order and the merged
   # tree's entry N no longer corresponded to raw event N. That silently breaks every positional
   # join: the IC join, and any (run, event, trackID) match to the fits. The PID plane is a
   # histogram and does not care, which is exactly why it went unnoticed.
   parts=()
   for ((c=0; c<NCHUNK; c++)); do
      f="$RECO/.chunk/$(printf c%02d_ "$c")${run}_reco.root"
      [[ -s "$f" ]] && parts+=("$f")
   done
   if (( fail )) || [[ ! -s "${parts[0]}" ]]; then
      # RETRY ONCE, SERIALLY. The failure mode above is a concurrent-open race, so a single
      # reader almost always succeeds where 12 did not. Slower, but it recovers the run instead
      # of losing it -- and only the runs that actually failed pay the cost.
      echo "[$run] chunked pass failed -- retrying serially (one reader)"
      rm -f "$RECO"/.chunk/c[0-9][0-9]_"${run}"_reco.root
      if root -b -q "$HERE/unpackReco_C15d.C(\"$run\",$nev,false,\"$RECO/.chunk/c00_\",\"$RAW\",\"$PAR\",\"$GEO\",20.0,15.0,7.5,kTRUE,kFALSE,\"\",30,0.0,2.85,0)" \
            2>&1 | head -c 50000000 >"$LOG/${run}_serial.log" \
         && [[ -s "$RECO/.chunk/c00_${run}_reco.root" ]]; then
         parts=("$RECO/.chunk/c00_${run}_reco.root")
         echo "[$run] serial retry OK"
      else
         echo "[$run] FAILED even serially (see $LOG/${run}_serial.log)"
         rm -f "$RECO"/.chunk/c[0-9][0-9]_"${run}"_reco.root
         echo "$run" >> "$HERE/runs_reco_failed.txt"
         continue
      fi
   fi
   if hadd -f "$RECO/.chunk/${run}_merged.root" "${parts[@]}" >"$LOG/${run}_hadd.log" 2>&1 \
      && [[ -s "$RECO/.chunk/${run}_merged.root" ]]; then
      mv "$RECO/.chunk/${run}_merged.root" "$out"
      rm -f "$RECO"/.chunk/c[0-9][0-9]_"${run}"_reco.root
      # ★ The PID NTUPLE is a separate step and it is NOT optional. mkpid_C15d.C,
      # make_points_C15d.C and gain_report_C15d.C all read <run>_pid.root, not the reco. Omitting
      # it here (reco_batch.sh does run it) left the gain report with "runs with data : 1" and no
      # usable PID plane, which looks like a data problem rather than a missing stage.
      root -b -q "$HERE/pidntuple_C15d.C(\"$run\",\"$RECO/\")" >>"$LOG/${run}_pid.log" 2>&1 || true
      echo "[$run] OK  $nev events, ${#parts[@]} chunks, $((SECONDS-t0)) s (prefetch ${t_pf} s)  $(du -h "$out" | cut -f1)"
   else
      echo "[$run] HADD FAILED (see $LOG/${run}_hadd.log)"
      rm -f "$RECO"/.chunk/c[0-9][0-9]_"${run}"_reco.root
   fi
done
echo "=== done: $(ls -1 "$RECO"/*_reco.root 2>/dev/null | wc -l) reco files ==="
