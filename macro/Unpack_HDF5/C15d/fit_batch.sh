#!/usr/bin/env bash
# C15d: GENFIT+CATIMA fitting over the D2 run set, one species hypothesis per pass. Resumable.
#
#   ./fit_batch.sh [nparallel] [species] [runlist] [nEvents]
#     species: d (deuteron, default) | p (proton) | t (triton)
#
# Per run: fit -> <run>_genfit_<sp>.root, then flatten to <run>_kin_<sp>.root. Both are skipped if
# they already exist, so re-running continues where it stopped.
#
# ★ THE FIT IS UNGATED AND THE KIN NTUPLE HOLDS EVERY FITTED TRACK. The PID gate is applied later,
# by joining the kin ntuple on (run, event, trackID). Two reasons:
#   - AtGenfitter::SetPIDGate runs its OWN AtSpyralPID on RAW dE/dx, while every gate in this
#     workspace is drawn on the gain-matched plane. An in-fit gate would silently select a
#     different set of tracks than the one that was drawn.
#   - fitting once and gating many times means a gate can be revised without refitting.
#
# ★ USE GENFIT+CATIMA, NOT THE UKF (project decision). Material effects are ON by default in
# fitGenfit_C15d.C with the CATIMA MSC, straggling and dE/dx backends, and matFallback OFF so a
# failed material-effects fit drops out rather than being silently refitted without them.
#
# Bz = -2.85 T, measured: at +2.85 the deuterons come out at 0.04 MeV with chi2/ndf 4.0 against
# 2.02 MeV and 0.076 at -2.85.
#
# Sizes: a fit is ~160 MB per full run (~12 GB for the set) and a kin ntuple ~1.5 MB, so the fits
# are kept -- re-dumping with different cuts then costs nothing.

set -eo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/../../.." && pwd)"

NPAR="${1:-8}"
SPECIES="${2:-d}"
RUNLIST="${3:-$HERE/runs_d2.txt}"
NEVENTS="${4:--1}"

case "$SPECIES" in
   p) PDG=2212;       MASS=1.00782503207; ZED=1;;
   d) PDG=1000010020; MASS=2.01410178;    ZED=1;;
   t) PDG=1000010030; MASS=3.01550072;    ZED=1;;
   *) echo "usage: $0 [nparallel] [d|p|t] [runlist] [nEvents]" >&2; exit 1;;
esac

RECO_DIR="${C15D_RECO:-/home/yassid/C15d_reco}"
FIT_DIR="${C15D_FIT:-/home/yassid/C15d_fit}"
LOG_DIR="${C15D_LOGS:-/home/yassid/C15d_logs}"
MIN_FREE_GB="${MIN_FREE_GB:-40}"
BFIELD="${BFIELD:--2.85}"

set +u
# shellcheck disable=SC1091
source "$REPO/build/config.sh" >/dev/null 2>&1
set -u

[[ -n "${VMCWORKDIR:-}" ]] || { echo "ERROR: config.sh did not set VMCWORKDIR" >&2; exit 1; }
[[ -f "$RUNLIST" ]] || { echo "ERROR: run list not found: $RUNLIST" >&2; exit 1; }
mkdir -p "$FIT_DIR" "$LOG_DIR" "$FIT_DIR/.part"

mapfile -t RUNS < <(grep -vE '^\s*(#|$)' "$RUNLIST")
[[ ${#RUNS[@]} -gt 0 ]] || { echo "ERROR: empty run list -- refusing to no-op silently" >&2; exit 1; }

echo "=== C15d fit_batch  (GENFIT + CATIMA) ==="
echo "  species  : $SPECIES  (pdg $PDG, m $MASS u, Z $ZED)"
echo "  runs     : ${#RUNS[@]} from $(basename "$RUNLIST")"
echo "  parallel : $NPAR"
echo "  B        : $BFIELD T"
echo "  out      : $FIT_DIR"
echo "  free     : $(df -BG --output=avail "$FIT_DIR" | tail -1 | tr -d ' G') GB (guard ${MIN_FREE_GB})"
echo

trap 'pkill -P $$ 2>/dev/null || true' EXIT INT TERM

do_fit() {
   local run="$1"
   local fit="$FIT_DIR/${run}_genfit_${SPECIES}.root"
   local kin="$FIT_DIR/${run}_kin_${SPECIES}.root"
   local log="$LOG_DIR/${run}_fit_${SPECIES}.log"
   local reco="$RECO_DIR/${run}_reco.root"

   [[ -s "$reco" ]] || { echo "[$run] no reco, skipping"; return 0; }

   local free_gb
   free_gb=$(df -BG --output=avail "$FIT_DIR" | tail -1 | tr -d ' G')
   if (( free_gb < MIN_FREE_GB )); then echo "[$run] SKIP: only ${free_gb} GB free"; return 0; fi

   if [[ -s "$fit" ]]; then
      echo "[$run] fit exists, skipping"
   else
      # Stage and move: ROOT will "recover" a truncated fit and report plausible kinematics, so a
      # killed job must never leave one under its final name.
      local part="$FIT_DIR/.part/$run"
      rm -rf "$part"; mkdir -p "$part"
      if root -b -q "$HERE/fitGenfit_C15d.C(\"$run\", $NEVENTS, \"$RECO_DIR/\", \"\", \"$part/\", $BFIELD, 2, 5, \"\", 4.0, 5.0, 178.0, kTRUE, kTRUE, $PDG, $MASS, $ZED, \"$SPECIES\")" \
            >"$log" 2>&1 && [[ -s "$part/${run}_genfit_${SPECIES}.root" ]]; then
         mv -f "$part/${run}_genfit_${SPECIES}.root" "$fit"
         rm -rf "$part"
         echo "[$run] fit OK  ($(du -h "$fit" | cut -f1))"
      else
         echo "[$run] FIT FAILED -- see $log"
         rm -rf "$part"
         return 0
      fi
   fi

   if [[ -s "$kin" ]]; then
      echo "[$run] kin exists, skipping"
   else
      # No gate here on purpose: every fitted track goes in, and gates are applied downstream.
      if root -b -q "$HERE/dump_kine_C15d.C(\"$run\",\"\",\"$FIT_DIR/\",\"\",\"$SPECIES\")" >>"$log" 2>&1 \
            && [[ -s "$kin" ]]; then
         echo "[$run] kin OK  ($(grep -oP '\d+(?= written)' "$log" | tail -1) tracks)"
      else
         echo "[$run] KIN DUMP FAILED -- see $log"
      fi
   fi
}
export -f do_fit
export HERE RECO_DIR FIT_DIR LOG_DIR MIN_FREE_GB NEVENTS SPECIES PDG MASS ZED BFIELD

printf '%s\n' "${RUNS[@]}" | xargs -P "$NPAR" -I{} bash -c 'do_fit "$@"' _ {}

echo
echo "=== done: $(ls -1 "$FIT_DIR"/*_genfit_${SPECIES}.root 2>/dev/null | wc -l) fits, $(ls -1 "$FIT_DIR"/*_kin_${SPECIES}.root 2>/dev/null | wc -l) kin ntuples ==="
rmdir "$FIT_DIR/.part" 2>/dev/null || true
