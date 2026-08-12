#!/usr/bin/env bash
# A/B the memory cost of run->SetStoreTraj in the generation step.
#
# The generation stage climbed to 4.6 GB in 3 minutes with no plateau while the reco stage sat
# flat at 1 GB, and SetStoreTraj(kTRUE) is the obvious candidate: it keeps every Geant4 step point
# of every track, deltas included, for the whole run. Nothing downstream of the generator reads
# trajectories -- digitisation takes AtMCPoint, PSA takes pads, PID takes hits -- so if this is
# the cost, it is being paid for nothing outside the event display.
#
# 2000 events rather than 12000: the growth is linear, so the slope shows up in a sixth of the
# time. Peak RSS is sampled at 1 s, which is fine against a ~30 s run.
#
#   ./storetraj_ab.sh          # runs A (as-is) then B (patched to kFALSE), restores the macro
set -uo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MACRO="$DIR/C16_pp_a1975_sim.C"
OUT=/tmp/claude-1000/-home-yassid/5b3ae710-1d32-4ced-a4d7-f9be2f7654b0/scratchpad
NEV=2000
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
mkdir -p "$OUT"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"

cp "$MACRO" "$OUT/C16_pp_a1975_sim.C.orig"
trap 'cp "$OUT/C16_pp_a1975_sim.C.orig" "$MACRO"; echo "macro restored"' EXIT

run_one() { # $1 = label
   local label=$1 peak=0 rss
   cd "$DIR"
   root -b -q -l "C16_pp_a1975_sim.C($NEV,2.,178.,\"TGeant4\",-28.5,\"$OUT/ab_${label}.root\",0.0,7777)" \
      > "$OUT/ab_${label}.log" 2>&1 &
   local rpid=$!
   # Sample the ROOT child, not the wrapper: `root` is a shell script that execs root.exe.
   while kill -0 $rpid 2>/dev/null; do
      for p in $(pgrep -P $rpid 2>/dev/null) $rpid; do
         rss=$(awk '/^VmRSS:/{print $2}' "/proc/$p/status" 2>/dev/null)
         [ -n "${rss:-}" ] && [ "$rss" -gt "$peak" ] && peak=$rss
      done
      sleep 1
   done
   wait $rpid; local rc=$?
   local secs=$(grep -oP 'Real time \K[0-9.]+' "$OUT/ab_${label}.log" | tail -1)
   printf '%-28s peak RSS %6d MB   real %s s   rc=%d\n' "$label" $((peak/1024)) "${secs:-?}" $rc
}

echo "=== A: SetStoreTraj as committed (kTRUE) ==="
run_one storetraj_true

echo "=== B: SetStoreTraj(kFALSE) ==="
sed -i 's/run->SetStoreTraj(kTRUE);/run->SetStoreTraj(kFALSE);/' "$MACRO"
grep -n "SetStoreTraj" "$MACRO"
run_one storetraj_false

echo "=== output sizes ==="
ls -l "$OUT"/ab_storetraj_*.root 2>/dev/null | awk '{print $9, $5/1024/1024" MB"}'
