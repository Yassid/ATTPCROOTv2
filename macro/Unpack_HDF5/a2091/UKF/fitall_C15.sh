#!/bin/bash
# ============================================================================
# Standalone Phase C for a2091 15C+p: ungated three-fitter pass over every
# reco'd run.  UKF (build/) + GENFIT matEffects=OFF + GENFIT matEffects=ON
# (both build_genfit/), all over the SAME reco, so the three are comparable.
#
# This replaces overnight_C15.sh's Phase C, which could never run: that script
# uses `set -u`, and sourcing thisroot.sh under `set -u` kills the shell
# (unbound variable, exit 127, stderr sent to /dev/null -> silent death).
# Here the sources are wrapped in `set +u`.
#
# Auto-discovers runs from $RECO, skips reco stubs (<5 MB: run_0140 corrupt,
# run_0180 junk, and any not-yet-reco'd FRIBDAQ run), and skips fits that
# already exist -- so it is safe to re-run to pick up newly reco'd runs.
#
#   ./fitall_C15.sh            # all reco'd runs, 6 parallel
#   ./fitall_C15.sh "run_0138 run_0179" 4
# ============================================================================
NPAR="${2:-6}"
REPO="/home/yassid/fair_install/ATTPCROOTv2"
HERE="$REPO/macro/Unpack_HDF5/a2091/UKF"
RECO="${RECO:-/home/yassid/a2091_C15_reco}"
FIT="${FIT:-/home/yassid/a2091_C15_fit}"
LOG="$FIT/logs"; mkdir -p "$FIT" "$LOG"
MASTER="$FIT/fitall.log"
MIN_FREE_GB="${MIN_FREE_GB:-15}"
GEOM="${GEOM:-ATTPC_H300torr_RT}"
DENS="${DENS:-3.308e-5}"          # H2, 300 torr, room temperature
MINRECO=5000000                   # bytes; below this the reco is a stub

# Killing this script must also take down its xargs pool -- otherwise xargs keeps
# spawning fits after the parent is gone (learned the hard way).
trap 'pkill -P $$ 2>/dev/null' EXIT INT TERM

say(){ echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$MASTER"; }
# Free space on the filesystem we actually WRITE to -- not always /home, since $FIT
# can point at the external drive. Hardcoding /home made the guard veto every fit
# when /home was full but the output disk had 1.7 TB free.
free_gb(){ df -BG --output=avail "$FIT" 2>/dev/null | tail -1 | tr -dc '0-9'; }
guard(){ local f; f=$(free_gb)
         if [ "${f:-0}" -lt "$MIN_FREE_GB" ]; then
           echo "DISK GUARD: ${f}GB free < ${MIN_FREE_GB}GB -- skipping $1"; return 1; fi; return 0; }

# Load a build's environment. thisroot.sh and config.sh both trip `set -u`,
# hence the explicit set +u; also both are noisy, hence the redirects.
srcbuild(){ set +u
            source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
            source "$REPO/$1/config.sh" >/dev/null 2>&1
            export ROOT_INCLUDE_PATH="$REPO/$1/include:/home/yassid/fair_install/FairRootInstall/include"; }

# A reco that is still being written must never be fit: it would be fit over a
# truncated tree and look like a low-statistics run. Two independent checks --
# no live unpackReco for that run, and the file untouched for MINAGE seconds.
MINAGE="${MINAGE:-600}"
is_stable(){ local f="$1" r="$2"
  # pgrep, not `ps | grep`: grep's own command line contains the pattern and so
  # matches itself, which made every run look like it was still being written.
  # Parens/dots escaped because pgrep -f takes an ERE.
  pgrep -f "unpackReco_C15\.C\(\"$r\"" >/dev/null 2>&1 && return 1
  local age=$(( $(date +%s) - $(stat -c%Y "$f") ))
  [ "$age" -ge "$MINAGE" ]
}

# Runs to fit: explicit list, else every complete, sufficiently-large reco in $RECO.
if [ -n "${1:-}" ]; then
  RUNS="$1"
else
  RUNS=$(for f in "$RECO"/*_reco.root; do
           [ -f "$f" ] || continue
           [ "$(stat -c%s "$f")" -lt "$MINRECO" ] && continue
           r=$(basename "$f" _reco.root)
           is_stable "$f" "$r" || { echo "SKIP $r (reco still being written)" >&2; continue; }
           echo "$r"
         done | sort -u | tr '\n' ' ')
fi
say "=== fitall 15C+p: $(echo $RUNS | wc -w) runs, NPAR=$NPAR, free=$(free_gb)GB ==="
say "    geom=$GEOM  density=$DENS  out=$FIT"

# ---- one run, one fitter -----------------------------------------------------
ukf_one(){ local r="$1"
  [ -f "$FIT/${r}_ukf.root" ] && { echo "skip $r (ukf exists)"; return; }
  guard "$r ukf" || return
  root -b -q -l "$HERE/pipeline/fitUKF_C15.C(\"$r\",-1,\"proton\",-1,2.85,$DENS,\"\",\"$RECO/\",0.5,0.1,1,10,\"$FIT/\")" \
       > "$LOG/${r}_ukf.log" 2>&1
  echo "[$(date +%H:%M:%S)] $r ukf=$([ -f "$FIT/${r}_ukf.root" ] && echo ok || echo FAIL)"
}
genfit_one(){ local r="$1" tag="$2" mat="$3"
  [ -f "$FIT/${r}_genfit${tag}.root" ] && { echo "skip $r (genfit$tag exists)"; return; }
  guard "$r genfit$tag" || return
  root -b -q -l "$HERE/pipeline/fitGenfit_C15.C(\"$r\",-1,\"$RECO/\",\"$tag\",\"$FIT/\",-2.85,2,5,\"\",4.0,10.0,170.0,$mat,kFALSE,\"proton\",\"$GEOM\")" \
       > "$LOG/${r}_genfit${tag}.log" 2>&1
  echo "[$(date +%H:%M:%S)] $r genfit$tag=$([ -f "$FIT/${r}_genfit${tag}.root" ] && echo ok || echo FAIL)"
}
export -f ukf_one genfit_one guard free_gb   # guard calls free_gb inside the xargs subshells
export HERE RECO FIT LOG DENS GEOM MIN_FREE_GB

# ---- C1: UKF, from build/ ----------------------------------------------------
srcbuild build
say "C1: UKF (build/)"
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'ukf_one "$@"' _ {} | tee -a "$MASTER"
say "C1 done: $(ls "$FIT"/*_ukf.root 2>/dev/null | wc -l) ukf files, free=$(free_gb)GB"

# ---- C2/C3: GENFIT, from build_genfit/ --------------------------------------
srcbuild build_genfit
say "C2: GENFIT matEffects=OFF (build_genfit/)"
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'genfit_one "$1" _nomat kFALSE' _ {} | tee -a "$MASTER"
say "C2 done: $(ls "$FIT"/*_genfit_nomat.root 2>/dev/null | wc -l) files, free=$(free_gb)GB"

say "C3: GENFIT matEffects=ON (build_genfit/)"
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'genfit_one "$1" _mat kTRUE' _ {} | tee -a "$MASTER"
say "C3 done: $(ls "$FIT"/*_genfit_mat.root 2>/dev/null | wc -l) files, free=$(free_gb)GB"

say "=== fitall DONE.  ukf=$(ls "$FIT"/*_ukf.root 2>/dev/null | wc -l)  nomat=$(ls "$FIT"/*_genfit_nomat.root 2>/dev/null | wc -l)  matFX=$(ls "$FIT"/*_genfit_mat.root 2>/dev/null | wc -l) ==="
