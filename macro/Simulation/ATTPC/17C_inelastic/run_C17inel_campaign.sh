#!/usr/bin/env bash
# The 17C M_n/M_p proposal simulation: three bound levels of 17C, two probes, two fields.
#
#   ./run_C17inel_campaign.sh [-j N] [nEvents]
#
# THE MATRIX. 2 channels x 3 levels x 2 fields = 12 samples.
#
#   channel   pp   17C(p,p') on 300 torr H2      the proposal's day 1
#             dd   17C(d,d') on 300 torr D2      the proposal's day 2
#   level     gs      0      3/2+                elastic
#             ex217   0.217  1/2+                proposal target
#             ex332   0.332  5/2+                proposal target
#   field     2.85 T  what Ref.[24] ran at, and this proposal's implicit configuration
#             4.00 T  SOLARIS's design field
#
# There is NO pad-pitch axis: this proposal runs the conventional AT-TPC pad plane, so the pitch
# comparison the 14C campaign made is not on offer here.
#
# WHY A FIELD AXIS AT ALL. inel_kinematics_C17.C settles it before any event is generated. The two
# proposal states are 115 keV apart, which appears as 217 keV of recoil energy at EVERY lab angle
# (dEx/dKE is flat at -0.53). sigma(Ex) is therefore set entirely by sigma(KE), which the 14C
# matrix MEASURED on this same pad plane as 0.343 MeV at 2.85 T and 0.067 MeV at 4 T. Propagated,
# that is sigma(Ex) = 0.18 MeV at 2.85 T (separation 0.32 -- one bump) against 0.043 MeV at 4 T
# (separation 1.34). The field is the only lever there is, and the campaign exists to confirm that
# prediction with a full simulation rather than a propagation.
#
# WAVES. Wave 1 is the whole (p,p') channel at both fields -- the primary probe, and by itself the
# complete answer to "can the doublet be decomposed". Wave 2 adds (d,d'), which the ratio needs but
# which cannot change the (p,p') conclusion. INEL_ONE_WAVE=1 merges them: on a resume whose
# generations are already on disk reco dominates, and two waves of six leave half the jobs idle.
set -uo pipefail
JOBS=8
if [ "${1:-}" = "-j" ]; then JOBS=${2:?-j needs a number}; shift 2; fi
NEV=${1:-16000}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOTDIR=${INEL_ROOT:-/media/yassid/Seagate Hub/ATTPC/C17_inel}
export INEL_ROOT="$ROOTDIR"
mkdir -p "$ROOTDIR"
MASTER="$ROOTDIR/campaign.log"
say() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$MASTER"; }

# ONE DRIVER PER OUTPUT TREE. Two campaigns sharing an INEL_ROOT do not merely duplicate work: they
# write the same reco and genfit files concurrently, and each one's resume check then sees the
# other's half-written output. That is how the (d,p) arm's first run failed, losing two hours. The
# stale-lock branch checks the recorded PID rather than just the directory, so a driver killed
# mid-run does not block the next one forever.
DRVLOCK="$ROOTDIR/campaign.lock"
if ! mkdir "$DRVLOCK" 2>/dev/null; then
   other=$(cat "$DRVLOCK/pid" 2>/dev/null || echo "?")
   if [ "$other" != "?" ] && kill -0 "$other" 2>/dev/null; then
      echo "REFUSING TO START: campaign driver PID $other is already running on $ROOTDIR"
      echo "  (kill it, or set INEL_ROOT to a different tree)"
      exit 2
   fi
   echo "removing stale lock from PID ${other:-?} (no such process)"
   rm -rf "$DRVLOCK"; mkdir "$DRVLOCK" || { echo "cannot take $DRVLOCK"; exit 2; }
fi
echo $$ > "$DRVLOCK/pid"
trap 'rm -rf "$DRVLOCK"' EXIT INT TERM

# channel:state:field:seed -- distinct seeds throughout, or parallel jobs generate byte-identical
# events and the added statistics are a copy of the same sample.
WAVE1="pp:gs:2.85:5101 pp:ex217:2.85:5102 pp:ex332:2.85:5103 \
       pp:gs:4.0:5104  pp:ex217:4.0:5105  pp:ex332:4.0:5106"
WAVE2="dd:gs:2.85:5201 dd:ex217:2.85:5202 dd:ex332:2.85:5203 \
       dd:gs:4.0:5204  dd:ex217:4.0:5205  dd:ex332:4.0:5206"

if [ "${INEL_ONE_WAVE:-0}" = "1" ]; then
   WAVE1="$WAVE1 $WAVE2"
   WAVE2=""
fi

run_wave() {
   local label=$1 list=$2
   say "=== wave $label : $(echo $list | wc -w) samples, -j $JOBS, $NEV events each ==="
   printf "%s\n" $list | xargs -P "$JOBS" -I{} bash -c '
     IFS=: read -r ch st bt seed <<< "{}"
     '"$HERE"'/accumulate_C17inel.sh "$ch" "$st" "$bt" "$seed" '"$NEV"' 2>&1 | tail -3'
   say "=== wave $label finished ==="
}

say "########## 17C(p,p')/(d,d') campaign start (-j $JOBS, $NEV events/sample, out $ROOTDIR) ##########"
run_wave 1 "$WAVE1"
[ -n "$WAVE2" ] && run_wave 2 "$WAVE2"
say "########## 17C(p,p')/(d,d') campaign done ##########"

n=$(ls "$ROOTDIR"/*/*.marker 2>/dev/null | wc -l)
say "completed samples: $n/12"
for cfg in pp_b285 pp_b400 dd_b285 dd_b400; do
   for st in gs ex217 ex332; do
      acc=$(grep -h 'overall acceptance' "$ROOTDIR/$cfg/"*_"${st}"_*_acc.log 2>/dev/null | tail -1)
      say "$(printf '%-10s %-6s %s' "$cfg" "$st" "${acc:-MISSING}")"
   done
done
