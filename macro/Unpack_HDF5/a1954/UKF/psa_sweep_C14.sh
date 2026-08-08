#!/usr/bin/env bash
# Does the +8 % proton KE excess above 8.5 MeV come from extra hits in the real point cloud?
#
# Space charge is excluded quantitatively (lambda ~ 6.5e-12 C/m at 2000 pps gives at most 0.23 mm
# of radial distortion, against a 0.5 mm hit sigma), so the remaining data-vs-sim asymmetry is the
# point cloud itself: the simulation carries ~1.8x fewer pads per track than the data. Extra hits
# -- noise, delta electrons, beam halo -- would pull a partial-arc fit exactly where it is weakest,
# which is above 8.5 MeV where the helix no longer closes inside the detector.
#
# The test varies only what enters the fit, on the SAME events:
#   thr 20 + cleaner  (baseline, reused from the triplclust/HDBSCAN test)
#   thr 40 + cleaner  (drop low-amplitude hits)
#   thr 80 + cleaner  (drop aggressively)
#   thr 20, NO cleaner (AtDirDeDxCleaner off -- the opposite direction)
# If the bias moves with the hit content, it is contamination. If it does not, the hits themselves
# are mispositioned and the next suspect is the PSA z-assignment.
#
#   ./psa_sweep_C14.sh run_0058 12000

set -uo pipefail
RUN="${1:-run_0058}"
NEV="${2:-12000}"
IN="/mnt/g/a1954_remerged/"
OUT="/home/yassid/a1954_C14_psasweep"
HERE="$(cd "$(dirname "$0")" && pwd)"
REPO="$(cd "$HERE/../../../.." && pwd)"
PAR="ATTPC.a1954_C14.par"

set +u
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
set -u

[[ -d "$IN" ]] || { echo "FATAL: $IN not found" >&2; exit 1; }

# variant | PSA threshold | doClean
VARIANTS=("thr40|40|kTRUE" "thr80|80|kTRUE" "noclean|20|kFALSE")

run_one() {
   local tag="$1" thr="$2" clean="$3"
   local d="$OUT/$tag"
   mkdir -p "$d/logs"
   if [[ -f "$d/${RUN}.marker" ]]; then echo "skip $tag (done)"; return 0; fi
   echo "=== $tag : thr=$thr doClean=$clean ==="
   root -b -q -l "$HERE/pipeline/unpackReco_C14.C(\"$RUN\",$NEV,false,\"$d/\",\"$IN\",false,false,\"multifit\",0,$thr,\"hdbscan\",20,8,0,0.1,\"$PAR\",$clean)" \
        > "$d/logs/${RUN}_reco.log" 2>&1
   [[ $? -eq 0 && -s "$d/${RUN}_reco.root" ]] || { echo "FAILED reco $tag" >&2; return 1; }
   root -b -q -l "$HERE/pipeline/fitUKF_C14.C(\"$RUN\",-1,\"proton\",-1,2.85,3.553e-5,\"\",\"$d/\",0.5,0.1,1,10,\"$d/\")" \
        > "$d/logs/${RUN}_fit.log" 2>&1
   [[ -s "$d/${RUN}_ukf.root" ]] || { echo "FAILED fit $tag" >&2; return 1; }
   root -b -q -l "$HERE/pp/ex_C14.C(\"$RUN\",\"$d/\",161.0,1e9,\"_sweep_${tag}\",1.007825,14.003242,\"\",\"ukf\",kFALSE)" \
        > "$d/logs/${RUN}_ex.log" 2>&1
   grep -E "good track" "$d/logs/${RUN}_ex.log" | tail -1
   touch "$d/${RUN}.marker"          # marker means COMPLETED, not "a file exists"
   echo "  done $tag"
}

# two at a time: each reco peaks around 2-3 GB and the box has ~20 GB free
pids=()
for v in "${VARIANTS[@]}"; do
   IFS='|' read -r tag thr clean <<< "$v"
   run_one "$tag" "$thr" "$clean" &
   pids+=($!)
   while [[ $(jobs -rp | wc -l) -ge 2 ]]; do sleep 20; done
done
wait "${pids[@]}"
echo "SWEEP RECO+FIT DONE"
