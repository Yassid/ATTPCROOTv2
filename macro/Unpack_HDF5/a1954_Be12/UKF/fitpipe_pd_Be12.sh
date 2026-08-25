#!/bin/bash
# 12Be(p,d)11Be pass: per run  gate (IC 500-900 + DEUTERON-gated tracks) -> UKF(deuteron).
# Separate dirs from the (p,p') pass so nothing is overwritten. GENFIT is skipped on
# purpose (UKF beats it on this data set, see ANALYSIS_REPORT).
#   ./fitpipe_pd_Be12.sh "run_0143 run_0147" 4 [thMin]
RUNS="${1:-run_0143}"; NPAR="${2:-4}"; THMIN="${3:-0}"; MP="${4:-30}"
# 12Be beam window on the FRIB IC amplitude: 500-800 (Yassid, 2026-08-25), replacing the 500-900
# this pipeline used in July. MP must match the plane the deuteron gate was drawn on.
ICLO="${ICLO:-500}"; ICHI="${ICHI:-800}"
REPO="/home/yassid/fair_install/ATTPCROOTv2-OpenKF"
HERE="$REPO/macro/Unpack_HDF5/a1954_Be12/UKF"
SLIM="/home/yassid/a1954_Be12_reco_hdb_slim/"
FREF="/mnt/h/a1954_Be12_reco_hdb/"   # the F: copy is gone
FITDIR="/home/yassid/a1954_Be12_fit_pd/"; IN="${FITDIR}in/"; LOG="${FITDIR}logs"; mkdir -p "$IN" "$LOG"
DGATE="$HERE/pid/deuteron_12Be.json"

source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"

# The gated input is reused ONLY if the configuration matches -- see fitpipe_Be12.sh.
gatecfg(){ printf 'ic=%s-%s tb=%s-%s th=%s mp=%s gate=%s slim=%s\n' \
             "$ICLO" "$ICHI" 1050 1250 "$THMIN" "$MP" "$(md5sum "$DGATE" 2>/dev/null | cut -d" " -f1)" "$SLIM"; }
export -f gatecfg

one(){ local r="$1"; local L="$LOG/${r}.log"; local stamp="$IN/${r}.gatecfg"; : > "$L"
  local want="$(gatecfg)"
  if [ -f "$IN/${r}_reco.root" ] && [ "$(cat "$stamp" 2>/dev/null)" = "$want" ]; then
    echo "gate: reusing $IN/${r}_reco.root  [$want]" >> "$L"
  else
    [ -f "$SLIM/${r}_FRIB.root" ] || { echo "no FRIB $r"; return; }
    [ $(stat -c%s "$SLIM/${r}_FRIB.root" 2>/dev/null) -lt 10000 ] && { echo "empty FRIB $r (skip)"; return; }
    if [ -f "$IN/${r}_reco.root" ]; then
      echo "gate: REBUILDING -- config changed" >> "$L"
      echo "    was: $(cat "$stamp" 2>/dev/null || echo '(no stamp: predates this check)')" >> "$L"
      echo "    now: $want" >> "$L"
      rm -f "$IN/${r}_reco.root"
    fi
    root -b -q -l "$HERE/pipeline/gate_events_Be12.C(\"$r\",\"$SLIM\",\"$IN\",\"$FREF${r}_reco.root\",$ICLO,$ICHI,1050,1250,2.85,\"$DGATE\",$THMIN,200,800,1500,$MP,0.0)" >> "$L" 2>&1
    [ -f "$IN/${r}_reco.root" ] && printf '%s' "$want" > "$stamp"
  fi
  [ -f "$IN/${r}_reco.root" ] || { echo "gate failed $r"; return; }
  # 2) UKF with the DEUTERON mass hypothesis
  root -b -q -l "$HERE/pipeline/fitUKF_Be12.C(\"$r\",-1,\"deuteron\",-1,2.85,3.308e-5,\"\",\"$IN\",0.5,0.1,1,10,\"$FITDIR\")" >> "$L" 2>&1
  # 3) GENFIT, deuteron, material effects + CATIMA (defaults), forward theta window (5 deg).
  #    July skipped this because the UKF beat it -- but that was matEffects=kFALSE against a
  #    2x-too-dense geometry. Both fitters now go into the SAME dir so the explorer can switch.
  root -b -q -l "$HERE/pipeline/fitGenfit_Be12.C(\"$r\",-1,\"$IN\",\"\",\"$FITDIR\",-2.85,2,5,\"\",4.0,5.0,170.0,kTRUE,kFALSE,\"deuteron\")" >> "$L" 2>&1
  echo "[$(date +%H:%M:%S)] $r  ukf=$([ -f $FITDIR${r}_ukf.root ]&&echo ok||echo FAIL)  genfit=$([ -f $FITDIR${r}_genfit.root ]&&echo ok||echo FAIL)"
}
export -f one; export HERE SLIM FREF FITDIR IN LOG DGATE THMIN MP ICLO ICHI
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "[$(date +%H:%M:%S)] (p,d) FIT PIPE DONE -> $FITDIR"
