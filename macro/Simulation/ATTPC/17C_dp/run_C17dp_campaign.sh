#!/usr/bin/env bash
# The 17C(d,p)18C proposal simulation: the four bound states of 18C at the nominal
# SOLARIS + AT-TPC configuration (D2 300 torr, B = 2.85 T, real AT-TPC pad plane).
#
#   ./run_C17dp_campaign.sh [-j N] [nEvents]
#
# Four samples, one per level. There is no field or pad-pitch axis: this is a proposal simulation
# of the detector as it exists, not the design matrix that 14C(d,p) was. The 14C(d,p) matrix
# already answered the design question for a backward-ejectile transfer at this beam energy --
# backward sigma(Ex) is FLAT at ~0.2 MeV across 2.85/4/7 T and both pad planes, so neither field
# nor pitch buys anything where the transfer peaks (see 14C_dp/RESULTS.md).
#
# 18C has S_n = 4184 keV, so these four levels ARE the bound spectrum:
#     gs      0        0+
#     ex1588  1.588    2+          the state the proposal's B(E2) discussion is about
#     ex2515  2.515    (2+)
#     ex3972  3.972    (2,3)+
# Wave 1 is the two lowest, which carry the 1588 keV and 927 keV separations that decide whether
# the spectrum is resolvable at all; wave 2 adds the upper pair.
set -uo pipefail
JOBS=4
if [ "${1:-}" = "-j" ]; then JOBS=${2:?-j needs a number}; shift 2; fi
NEV=${1:-12000}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
ROOTDIR=${DP_ROOT:-/mnt/f/C17dp}
export DP_ROOT="$ROOTDIR"
mkdir -p "$ROOTDIR"
MASTER="$ROOTDIR/campaign.log"
say() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$MASTER"; }

# ONE DRIVER PER OUTPUT TREE. Two campaigns sharing a DP_ROOT do not merely duplicate work: they
# write the same reco and genfit files concurrently, and each one's resume check then sees the
# other's half-written output. That is exactly how the first run of this campaign failed -- a
# driver whose children had been killed was still alive, a second was started, and the second read
# a genfit file the first was still writing ("probably not closed / missing cbmsim"), losing two
# hours. The per-sample generation lock did not help, because the collision was downstream of it.
# The stale-lock branch checks the recorded PID rather than just the directory, so a driver killed
# mid-run does not block the next one forever.
DRVLOCK="$ROOTDIR/campaign.lock"
if ! mkdir "$DRVLOCK" 2>/dev/null; then
   other=$(cat "$DRVLOCK/pid" 2>/dev/null || echo "?")
   if [ "$other" != "?" ] && kill -0 "$other" 2>/dev/null; then
      echo "REFUSING TO START: campaign driver PID $other is already running on $ROOTDIR"
      echo "  (kill it, or set DP_ROOT to a different tree)"
      exit 2
   fi
   echo "removing stale lock from PID ${other:-?} (no such process)"
   rm -rf "$DRVLOCK"; mkdir "$DRVLOCK" || { echo "cannot take $DRVLOCK"; exit 2; }
fi
echo $$ > "$DRVLOCK/pid"
trap 'rm -rf "$DRVLOCK"' EXIT INT TERM

# state:seed -- distinct seeds, or parallel jobs generate byte-identical events and the added
# statistics are a copy of the same sample.
WAVE1="gs:9001 ex1588:9002"
WAVE2="ex2515:9003 ex3972:9004"

# The two-wave split exists so the two lowest levels -- which carry the 1588 keV and 927 keV
# separations that decide whether the spectrum is resolvable at all -- produce a headline number
# before the upper pair is finished. That is worth it on a cold start, where generation dominates.
# On a RESUME whose generations are already on disk, reco dominates instead, and splitting four
# samples into two waves of two leaves half the requested jobs idle. DP_ONE_WAVE=1 merges them.
if [ "${DP_ONE_WAVE:-0}" = "1" ]; then
   WAVE1="$WAVE1 $WAVE2"
   WAVE2=""
fi

run_wave() {
   local label=$1 list=$2
   say "=== wave $label : $(echo $list | wc -w) samples, -j $JOBS, $NEV events each ==="
   printf "%s\n" $list | xargs -P "$JOBS" -I{} bash -c '
     IFS=: read -r st seed <<< "{}"
     '"$HERE"'/accumulate_C17dp.sh "$st" "$seed" '"$NEV"' 2>&1 | tail -3'
   say "=== wave $label finished ==="
}

say "########## 17C(d,p)18C campaign start (-j $JOBS, $NEV events/sample, out $ROOTDIR) ##########"
run_wave 1 "$WAVE1"
[ -n "$WAVE2" ] && run_wave 2 "$WAVE2"
say "########## 17C(d,p)18C campaign done ##########"

n=$(ls "$ROOTDIR/b285_attpc"/*.marker 2>/dev/null | wc -l)
say "completed samples: $n/4"
for st in gs ex1588 ex2515 ex3972; do
   acc=$(grep -h 'overall acceptance' "$ROOTDIR/b285_attpc/${st}"_*_acc.log 2>/dev/null | tail -1)
   say "$(printf '%-8s %s' "$st" "${acc:-MISSING}")"
done
