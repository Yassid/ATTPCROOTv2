#!/usr/bin/env bash
# C15d: GENFIT+CATIMA fitting over the D2 run set, one species hypothesis per pass. Resumable.
#
#   ./fit_batch.sh [nparallel] [species] [runlist] [nEvents] [gate] [backExtrap]
#     species: d (deuteron, default) | p (proton) | t (triton)
#     gate   : a PID gate JSON -> the run is reduced BEFORE fitting. Empty = fit everything.
#
# Per run: fit -> <run>_genfit_<sp>.root, then flatten to <run>_kin_<sp>.root. Both are skipped if
# they already exist, so re-running continues where it stopped.
#
# ★ GATE BEFORE FITTING. With a gate, pid/gate_events_C15d.C first reduces the run to IC-passing
# events holding only gated tracks, and the fit runs on that: 617 proton tracks on run_0026 against
# 17,043 ungated, i.e. ~28x less fitting. The (d,d') pass was run ungated and did not need to be.
#
# The gate is NOT passed to AtGenfitter::SetPIDGate, which runs its own AtSpyralPID on RAW dE/dx
# while every gate here is drawn on the gain-matched plane -- measured, that selects 4,217 tracks
# where the plane selects 2,606. gate_events_C15d.C instead tests the polygon against the persisted
# AtPIDEvent, which IS the plane, and reproduces it to the track.
#
# The kin ntuple still holds every track that was fitted, so downstream cuts stay free.
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
# Gate BEFORE fitting. With a gate, each run is first reduced to IC-passing events holding only
# gated tracks (pid/gate_events_C15d.C), and the fit runs on that -- roughly an order of magnitude
# less fitting. Empty = fit everything, which is what the (d,d') pass did.
GATE="${5:-}"
BACKEXTRAP="${6:-kTRUE}"

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
echo "  gate     : ${GATE:-none (fitting everything)}"
echo "  backExtr : $BACKEXTRAP"
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
   local recoDir="$RECO_DIR/"

   [[ -s "$reco" ]] || { echo "[$run] no reco, skipping"; return 0; }

   # ---- optional pre-fit gating -------------------------------------------------------------
   if [[ -n "$GATE" ]]; then
      local gdir="$FIT_DIR/in_${SPECIES}"
      local gated="$gdir/${run}_reco.root"
      if [[ ! -s "$gated" ]]; then
         mkdir -p "$gdir"
         root -b -q "$HERE/pid/gate_events_C15d.C(\"$run\",\"$GATE\",\"$RECO_DIR/\",\"$gdir/\")" \
            >"$LOG_DIR/${run}_gate_${SPECIES}.log" 2>&1 || true
      fi
      if [[ ! -s "$gated" ]]; then
         echo "[$run] no gated events (see $LOG_DIR/${run}_gate_${SPECIES}.log)"
         return 0
      fi
      recoDir="$gdir/"
   fi

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
      if root -b -q "$HERE/fitGenfit_C15d.C(\"$run\", $NEVENTS, \"$recoDir\", \"\", \"$part/\", $BFIELD, 2, 5, \"\", 4.0, 5.0, 178.0, kTRUE, kTRUE, $PDG, $MASS, $ZED, \"$SPECIES\", \"_reco\", \"ATTPC_D300torr_v2_geomanager.root\", \"ATTPC.C15d_D2_300torr.par\", 6.5643e-5, 2, kFALSE, kTRUE, kTRUE, kTRUE, kFALSE, $BACKEXTRAP)" \
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
export HERE RECO_DIR FIT_DIR LOG_DIR MIN_FREE_GB NEVENTS SPECIES PDG MASS ZED BFIELD GATE BACKEXTRAP

printf '%s\n' "${RUNS[@]}" | xargs -P "$NPAR" -I{} bash -c 'do_fit "$@"' _ {}

echo
echo "=== done: $(ls -1 "$FIT_DIR"/*_genfit_${SPECIES}.root 2>/dev/null | wc -l) fits, $(ls -1 "$FIT_DIR"/*_kin_${SPECIES}.root 2>/dev/null | wc -l) kin ntuples ==="
rmdir "$FIT_DIR/.part" 2>/dev/null || true
