#!/bin/bash
# 12Be(p,t)10Be pass: per run  gate (IC + TRITON-gated tracks) -> UKF(triton) -> GENFIT(triton).
# Separate dirs from (p,p') and (p,d) so nothing is overwritten.
#
# WHY THIS CANNOT REUSE THE (p,p') GATED INPUTS: gate_events_Be12.C keeps only the tracks inside
# the PID polygon it is given, so a1954_Be12_fit/in/ holds PROTONS ONLY -- every triton was
# discarded at that step. (p,t) must be re-gated from the full reco.
#
# MP MUST MATCH THE PLANE THE GATE WAS DRAWN ON. AtSpyralPID::fMinPoints defaults to 30 and a gate
# drawn on an mp15 plane applied at mp30 is a DIFFERENT cut, silently.
#
#   ./fitpipe_pt_Be12.sh "run_0143 run_0147" 4 [thMin] [minPoints]
RUNS="${1:-run_0143}"; NPAR="${2:-4}"; THMIN="${3:-0}"; MP="${4:-30}"
# 12Be beam window on the FRIB IC amplitude: 500-800 (Yassid, 2026-08-25).
# This is NOT the 625-750 documented in pid/PID_COMPARISON.md, and NOT the 500-900 that the
# (p,p') and (p,d) drivers pass. The contaminant beam near IC 1900 outnumbers 12Be ~20:1 inside
# the triton gate, so this window is doing essentially all of the beam selection.
ICLO="${ICLO:-500}"; ICHI="${ICHI:-800}"
REPO="/home/yassid/fair_install/ATTPCROOTv2-OpenKF"
HERE="$REPO/macro/Unpack_HDF5/a1954_Be12/UKF"
SLIM="/home/yassid/a1954_Be12_reco_hdb_slim/"
FREF="/mnt/h/a1954_Be12_reco_hdb/"           # reco moved off F: ; the F: copy is gone
FITDIR="${PTDIR:-/home/yassid/a1954_Be12_fit_pt/}"; IN="${FITDIR}in/"; LOG="${FITDIR}logs"
mkdir -p "$IN" "$LOG"
TGATE="${TGATE:-$HERE/pid/triton_12Be.json}"
[ -f "$TGATE" ] || { echo "ERROR: no triton gate at $TGATE -- draw it first with pid/draw_gate_Be12.C"; exit 1; }

source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"

# Same rule as fitpipe_Be12.sh: the gated input is reused only if the CONFIGURATION matches.
gatecfg(){ printf 'ic=%s-%s tb=%s-%s th=%s mp=%s gate=%s slim=%s\n' \
             "$ICLO" "$ICHI" 1050 1250 "$THMIN" "$MP" "$(md5sum "$TGATE" 2>/dev/null | cut -d" " -f1)" "$SLIM"; }
export -f gatecfg

one(){ local r="$1"; local L="$LOG/${r}.log"
  local stamp="$IN/${r}.gatecfg"
  local want="$(gatecfg)"
  if [ -f "$IN/${r}_reco.root" ] && [ "$(cat "$stamp" 2>/dev/null)" = "$want" ]; then
    echo "gate: reusing $IN/${r}_reco.root  [$want]" > "$L"
  else
    [ -f "$SLIM/${r}_FRIB.root" ] || { echo "no FRIB $r"; return; }
    [ $(stat -c%s "$SLIM/${r}_FRIB.root" 2>/dev/null) -lt 10000 ] && { echo "empty FRIB $r (skip)"; return; }
    if [ -f "$IN/${r}_reco.root" ]; then
      echo "gate: REBUILDING -- config changed" > "$L"
      echo "    was: $(cat "$stamp" 2>/dev/null || echo '(no stamp)')" >> "$L"
      echo "    now: $want" >> "$L"
      rm -f "$IN/${r}_reco.root"
    fi
    # 1) gate: IC window + triton PID polygon, theta_lab > THMIN (tritons are FORWARD, so 0)
    root -b -q -l "$HERE/pipeline/gate_events_Be12.C(\"$r\",\"$SLIM\",\"$IN\",\"$FREF${r}_reco.root\",$ICLO,$ICHI,1050,1250,2.85,\"$TGATE\",$THMIN,200,800,1500,$MP,0.0)" >> "$L" 2>&1
    [ -f "$IN/${r}_reco.root" ] && printf '%s' "$want" > "$stamp"
  fi
  [ -f "$IN/${r}_reco.root" ] || { echo "gate failed $r"; return; }
  # 2) UKF with the TRITON mass hypothesis
  root -b -q -l "$HERE/pipeline/fitUKF_Be12.C(\"$r\",-1,\"triton\",-1,2.85,3.308e-5,\"\",\"$IN\",0.5,0.1,1,10,\"$FITDIR\")" >> "$L" 2>&1
  # 3) GENFIT, triton, material effects + CATIMA dE/dx (defaults), forward theta window
  root -b -q -l "$HERE/pipeline/fitGenfit_Be12.C(\"$r\",-1,\"$IN\",\"\",\"$FITDIR\",-2.85,2,5,\"\",4.0,3.0,170.0,kTRUE,kFALSE,\"triton\")" >> "$L" 2>&1
  echo "[$(date +%H:%M:%S)] $r  ukf=$([ -f $FITDIR${r}_ukf.root ]&&echo ok||echo FAIL)  genfit=$([ -f $FITDIR${r}_genfit.root ]&&echo ok||echo FAIL)"
}
export -f one; export HERE SLIM FREF FITDIR IN LOG TGATE THMIN MP ICLO ICHI
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "[$(date +%H:%M:%S)] (p,t) FIT PIPE DONE -> $FITDIR"
