#!/bin/bash
# Full-statistics Spyral PID plane for a2091 15C, built PER RUN IN PARALLEL.
#
# pid_C15.C over a 28-run CSV list is single-threaded and CPU-bound on the Spyral
# estimator (~10 min per big run -> 4-5 h serial). One process per run instead, then
# merge: pid_C15.C persists only the `spid` TNtuple, so the histograms are rebuilt
# from the merged ntuple by mkpid_C15.C.
#
#   ./pidall_C15.sh                 # every reco'd run, 8 parallel
#   ./pidall_C15.sh "run_0138 run_0179" 4
set -u
REPO=/home/yassid/fair_install/ATTPCROOTv2
HERE=$REPO/macro/Unpack_HDF5/a2091/UKF
RECO=/home/yassid/a2091_C15_reco
PLOTS=$HERE/pid/plots
PART=$PLOTS/parts
LOG=$RECO/logs/pid
NPAR="${2:-8}"
MINRECO=5000000
mkdir -p "$PART" "$LOG"

set +u
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
set -u

if [ -n "${1:-}" ]; then RUNS="$1"; else
  RUNS=$(for f in "$RECO"/*_reco.root; do
           [ -f "$f" ] || continue
           [ "$(stat -c%s "$f")" -ge "$MINRECO" ] && basename "$f" _reco.root
         done | sort -u | tr '\n' ' ')
fi
echo "[$(date +%H:%M:%S)] PID over $(echo $RUNS | wc -w) runs, NPAR=$NPAR"

one(){ local r="$1"
  # pid_C15.C writes $PLOTS/pid_C15<tag>.root; give each run its own tag, then relocate
  [ -f "$PART/${r}.root" ] && { echo "skip $r"; return; }
  root -b -q -l "$HERE/pid/pid_C15.C(\"$r\",\"$RECO/\",-1,2.85,\"_p_$r\")" > "$LOG/${r}.log" 2>&1
  if [ -f "$PLOTS/pid_C15_p_${r}.root" ]; then
    mv "$PLOTS/pid_C15_p_${r}.root" "$PART/${r}.root"
    rm -f "$PLOTS/pid_C15_p_${r}.png"
    echo "[$(date +%H:%M:%S)] $r ok ($(grep -o 'Spyral-valid=[0-9]*' "$LOG/${r}.log" | head -1))"
  else
    echo "[$(date +%H:%M:%S)] $r FAILED"
  fi
}
export -f one; export HERE RECO PLOTS PART LOG
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}

echo "[$(date +%H:%M:%S)] merging $(ls "$PART"/*.root 2>/dev/null | wc -l) parts"
hadd -f -k "$PLOTS/pid_C15_all.root" "$PART"/*.root > "$LOG/hadd.log" 2>&1 \
  && echo "[$(date +%H:%M:%S)] merged -> $PLOTS/pid_C15_all.root" \
  || { echo "HADD FAILED, see $LOG/hadd.log"; exit 1; }

root -b -q -l "$HERE/pid/mkpid_C15.C(\"$PLOTS/pid_C15_all.root\")" 2>&1 | tail -20
echo "[$(date +%H:%M:%S)] PID BUILD DONE"
