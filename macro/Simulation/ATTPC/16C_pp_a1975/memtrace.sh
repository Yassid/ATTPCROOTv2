#!/usr/bin/env bash
# Memory tracer + OOM guard for the PID accumulation.
#
# Two OOM shutdowns of the whole WSL VM happened during this job's reco stage, and both left no
# evidence behind: the VM dies, dmesg goes with it, and the only trace is a truncated ROOT file.
# So this samples memory to a file on F: (survives the VM) and, before the host gets into trouble,
# kills the ROOT child instead. accumulate_pid.sh is resumable per stage, so a killed stage costs
# that stage; a killed VM costs the session.
#
#   ./memtrace.sh <traceCsv> <stageFile> [guardPct] [periodSec]
#
# The guard only ever kills a root.exe whose command line names this job's directories -- an
# unrelated ROOT session of the user's is not this script's business.
set -uo pipefail
TRACE=${1:?need a trace file}; STAGEFILE=${2:?need a stage file}
GUARD=${3:-85}; PERIOD=${4:-5}
GUARDLOG="$(dirname "$TRACE")/memtrace_guard.log"
# Which ROOT processes belong to the job being guarded. Override with MINE=... to guard a
# different job -- the data refit (fitGenfitter_a1975) is a separate job that this script was
# not watching, and that is how the third VM shutdown happened.
MINE=${MINE:-'a1975_C16_pp_pid|C16_pp_a1975|pidPass_a1975'}

total_kb=$(awk '/^MemTotal:/{print $2}' /proc/meminfo)
[ -s "$TRACE" ] || echo "time,stage,used_mb,avail_mb,swap_used_mb,pct_used,njobs,sum_rss_mb,max_rss_mb,max_vsz_mb,max_pid,top_proc,top_rss_mb" > "$TRACE"

while :; do
   # MemAvailable, not free: page cache is reclaimable and counting it as used would fire the
   # guard on a job that is merely reading a 225 MB file.
   read -r avail_kb swaptot_kb swapfree_kb < <(awk '
      /^MemAvailable:/{a=$2} /^SwapTotal:/{st=$2} /^SwapFree:/{sf=$2}
      END{print a, st, sf}' /proc/meminfo)
   used_kb=$((total_kb - avail_kb))
   swap_kb=$((swaptot_kb - swapfree_kb))
   pct=$(( used_kb * 100 / total_kb ))

   # The job's ROOT processes. With several running, the sum is what threatens the host and the
   # largest is what the guard would kill, so both are recorded -- plus the count, because a
   # sample where a job has already died reads very differently from one where all are alive.
   procs=$(pgrep -f 'root\.exe' 2>/dev/null | while read -r p; do
            cmd=$(tr '\0' ' ' < "/proc/$p/cmdline" 2>/dev/null) || continue
            [[ "$cmd" =~ $MINE ]] || continue
            rss=$(awk '/^VmRSS:/{print $2}' "/proc/$p/status" 2>/dev/null); rss=${rss:-0}
            vsz=$(awk '/^VmSize:/{print $2}' "/proc/$p/status" 2>/dev/null); vsz=${vsz:-0}
            echo "$rss $vsz $p"
         done | sort -rn)
   njobs=$(printf '%s' "$procs" | grep -c . )
   sum_rss=$(printf '%s\n' "$procs" | awk '{s+=$1} END{print int(s/1024)}')
   sum_rss=${sum_rss:-0}
   job=$(printf '%s\n' "$procs" | head -1)
   if [ -n "$job" ]; then
      job_rss=$(( ${job%% *} / 1024 )); rest=${job#* }
      job_vsz=$(( ${rest%% *} / 1024 )); job_pid=${rest##* }
   else
      job_rss=0; job_vsz=0; job_pid=-
   fi

   top=$(ps -eo comm=,rss= --sort=-rss | head -1)
   stage=$(cat "$STAGEFILE" 2>/dev/null || echo "-")

   printf '%s,%s,%d,%d,%d,%d,%d,%d,%d,%d,%s,%s,%d\n' \
      "$(date +%H:%M:%S)" "$stage" $((used_kb/1024)) $((avail_kb/1024)) $((swap_kb/1024)) \
      "$pct" "$njobs" "$sum_rss" "$job_rss" "$job_vsz" "$job_pid" "${top%% *}" \
      $(( ${top##* } / 1024 )) >> "$TRACE"

   if [ "$pct" -ge "$GUARD" ] && [ "$job_pid" != "-" ]; then
      {
         echo "[$(date +%F' '%H:%M:%S)] GUARD FIRED at ${pct}% (>= ${GUARD}%), stage=$stage"
         echo "   killing pid $job_pid  rss=${job_rss} MB  vsz=${job_vsz} MB"
         tr '\0' ' ' < "/proc/$job_pid/cmdline" 2>/dev/null; echo
      } >> "$GUARDLOG"
      kill -9 "$job_pid" 2>/dev/null
   fi
   sleep "$PERIOD"
done
