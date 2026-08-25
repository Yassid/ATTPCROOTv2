#!/bin/bash
# Per run: gate (IC 500-900 events + proton-gated tracks) -> UKF fit -> GENFIT fit.
# 4 cores. Outputs local to FITDIR, then rsync to F. Input slim+FRIB local; metadata from F reco.
#   ./fitpipe_Be12.sh "run_0142 run_0143 ..." 4
RUNS="${1:-run_0142}"; NPAR="${2:-4}"
REPO="/home/yassid/fair_install/ATTPCROOTv2-OpenKF"
HERE="$REPO/macro/Unpack_HDF5/a1954_Be12/UKF"
HERE_PID="$HERE/pid"
# 12Be beam window on the FRIB IC amplitude: 500-800 (Yassid, 2026-08-25). The July production
# and this morning's 300-torr refit both used 500-900, which reaches past the 12Be peak into the
# structure near IC 1000. pid/PID_COMPARISON.md still documents 625-750; 500-800 supersedes both.
ICLO="${ICLO:-500}"; ICHI="${ICHI:-800}"
PGATE="${PGATE:-$HERE_PID/proton_12Be.json}"
THMIN="${THMIN:-90}"        # (p,p') protons are BACKWARD
MP="${MP:-30}"              # must match the plane any gate was drawn on
SLIM="/home/yassid/a1954_Be12_reco_hdb_slim/"
FREF="/mnt/h/a1954_Be12_reco_hdb/"   # moved off F: ; F: copy is gone
FITDIR="/home/yassid/a1954_Be12_fit/"; IN="${FITDIR}in/"; LOG="${FITDIR}logs"; mkdir -p "$IN" "$LOG"

source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"

# What the gated input depends on. If ANY of it changes the input must be rebuilt -- comparing
# only "does the file exist" is what let an IC-window change silently do nothing.
gatecfg(){ printf 'ic=%s-%s tb=%s-%s th=%s mp=%s gate=%s slim=%s\n' \
             "$ICLO" "$ICHI" 1050 1250 "$THMIN" "$MP" "$(md5sum "$PGATE" 2>/dev/null | cut -d" " -f1)" "$SLIM"; }
export -f gatecfg

one(){ local r="$1"; local L="$LOG/${r}.log"; local stamp="$IN/${r}.gatecfg"
  : > "$L"   # truncate: an appended-to stale log made a failed run look like a previous success
  if [ ! -f "$IN/${r}_reco.root" ]; then
    [ -f "$SLIM/${r}_FRIB.root" ] || { echo "no FRIB $r"; return; }
    [ $(stat -c%s "$SLIM/${r}_FRIB.root" 2>/dev/null) -lt 10000 ] && { echo "empty FRIB $r (skip)"; return; }
  fi
  # 1) gate -- reused ONLY if it was produced with the SAME configuration
  local want="$(gatecfg)"
  if [ -f "$IN/${r}_reco.root" ] && [ "$(cat "$stamp" 2>/dev/null)" = "$want" ]; then
    echo "[$(date +%H:%M:%S)] $r  gate: reusing $IN/${r}_reco.root  [$want]" > "$L"
  else
    if [ -f "$IN/${r}_reco.root" ]; then
      echo "[$(date +%H:%M:%S)] $r  gate: REBUILDING -- config changed" > "$L"
      echo "    was: $(cat "$stamp" 2>/dev/null || echo '(no stamp: predates this check)')" >> "$L"
      echo "    now: $want" >> "$L"
      rm -f "$IN/${r}_reco.root"
    fi
    root -b -q -l "$HERE/pipeline/gate_events_Be12.C(\"$r\",\"$SLIM\",\"$IN\",\"$FREF${r}_reco.root\",$ICLO,$ICHI,1050,1250,2.85,\"$PGATE\",$THMIN,200,800,1500,$MP,0.0)" >> "$L" 2>&1
    [ -f "$IN/${r}_reco.root" ] && printf '%s' "$want" > "$stamp"
  fi
  [ -f "$IN/${r}_reco.root" ] || { echo "gate failed $r"; return; }
  # 2) UKF
  root -b -q -l "$HERE/pipeline/fitUKF_Be12.C(\"$r\",-1,\"proton\",-1,2.85,3.308e-5,\"\",\"$IN\",0.5,0.1,1,10,\"$FITDIR\")" >> "$L" 2>&1
  # 3) GENFIT
  root -b -q -l "$HERE/pipeline/fitGenfit_Be12.C(\"$r\",-1,\"$IN\",\"\",\"$FITDIR\",-2.85,2,5,\"\",4.0,10.0,170.0,kTRUE,kFALSE,\"proton\")" >> "$L" 2>&1
  echo "[$(date +%H:%M:%S)] $r  gated=$(grep -o 'gated-proton events[^,]*' $L|head -1)  ukf=$([ -f $FITDIR${r}_ukf.root ]&&echo ok)  genfit=$([ -f $FITDIR${r}_genfit.root ]&&echo ok)"
}
export -f one; export HERE HERE_PID SLIM FREF FITDIR IN LOG ICLO ICHI PGATE THMIN MP
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}
echo "[$(date +%H:%M:%S)] FIT PIPE DONE. Syncing to F..."
mkdir -p /mnt/h/a1954_Be12_fit
rsync -a --include='*_ukf.root' --include='*_genfit.root' --exclude='*' "$FITDIR" /mnt/h/a1954_Be12_fit/ 2>/dev/null
echo "[$(date +%H:%M:%S)] SYNC DONE -> /mnt/h/a1954_Be12_fit/"
