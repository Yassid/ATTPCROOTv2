#!/bin/bash
# Gated production for a2091 15C(p,p'): IC(15C) + proton PID gate -> UKF fit -> Ex spectrum.
#
#   stage 1  gate_events_C15.C : IC window from pid/ic_15C.json (the 15C beam peak) AND the
#            proton gate pid/proton_15C.json, per track, plus polar > thMin. Writes a small
#            gated reco holding only the surviving proton tracks.
#   stage 2  fitUKF_C15.C on the gated reco.
#   stage 3  ex_C15.C at the CALIBRATED Ebeam (157 MeV, not the old 195 placeholder).
#
# Stage 1 is the expensive one: AtSpyralPID has to be recomputed per track because the runs
# were reco'd before AtPIDTask was added to the chain, so the PID is not stored in the files.
# Everything lands under the Seagate-symlinked fit directory.
set -u
REPO=/home/yassid/fair_install/ATTPCROOTv2
HERE=$REPO/macro/Unpack_HDF5/a2091/UKF
RECO=/home/yassid/a2091_C15_reco          # symlink -> Seagate
IC=/home/yassid/a2091_C15_ic              # symlink -> Seagate
GIN=/home/yassid/a2091_C15_fit/in         # gated reco  (symlink -> Seagate)
GFIT=/home/yassid/a2091_C15_fit/gated     # gated fits  (symlink -> Seagate)
LOG=/home/yassid/a2091_C15_fit/logs/gated
NPAR="${2:-6}"
EBEAM="${EBEAM:-157}"                     # MeV, from the elastic-ridge / tilt calibration
DENS="${DENS:-3.308e-5}"
MINRECO=5000000
mkdir -p "$GIN" "$GFIT" "$LOG"

trap 'pkill -P $$ 2>/dev/null' EXIT INT TERM

set +u
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
set -u

if [ -n "${1:-}" ]; then RUNS="$1"; else
  RUNS=$(for f in "$RECO"/*_reco.root; do
           [ -f "$f" ] || continue
           [ "$(stat -c%s "$f")" -ge "$MINRECO" ] || continue
           r=$(basename "$f" _reco.root)
           [ -f "$IC/${r}_ic.root" ] && echo "$r"     # need the IC summary to gate
         done | sort -u | tr '\n' ' ')
fi
echo "[$(date +%H:%M:%S)] gated pass over $(echo $RUNS | wc -w) runs, NPAR=$NPAR, Ebeam=$EBEAM MeV"

one(){ local r="$1"
  if [ ! -f "$GIN/${r}_reco.root" ]; then
    root -b -q -l "$HERE/pipeline/gate_events_C15.C(\"$r\",\"$RECO/\",\"$GIN/\",\"$RECO/${r}_reco.root\")" \
         > "$LOG/${r}_gate.log" 2>&1
  fi
  [ -f "$GIN/${r}_reco.root" ] || { echo "[$(date +%H:%M:%S)] $r GATE FAILED"; return; }
  if [ ! -f "$GFIT/${r}_ukf.root" ]; then
    root -b -q -l "$HERE/pipeline/fitUKF_C15.C(\"$r\",-1,\"proton\",-1,2.85,$DENS,\"\",\"$GIN/\",0.5,0.1,1,10,\"$GFIT/\")" \
         > "$LOG/${r}_ukf.log" 2>&1
  fi
  echo "[$(date +%H:%M:%S)] $r $(grep -o 'gated-proton events[^-]*' "$LOG/${r}_gate.log" 2>/dev/null | head -1) ukf=$([ -f "$GFIT/${r}_ukf.root" ] && echo ok || echo FAIL)"
}
export -f one; export HERE RECO IC GIN GFIT LOG DENS
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}

echo "[$(date +%H:%M:%S)] gated: $(ls "$GIN"/*_reco.root 2>/dev/null | wc -l) gated recos, $(ls "$GFIT"/*_ukf.root 2>/dev/null | wc -l) ukf fits"

# ---- stage 3: Ex spectrum at the calibrated beam energy ----
CSV=$(ls "$GFIT"/*_ukf.root 2>/dev/null | sed 's#.*/##;s/_ukf.root//' | sort -u | paste -sd,)
if [ -n "$CSV" ]; then
  echo "[$(date +%H:%M:%S)] Ex spectrum, gated, Ebeam=$EBEAM"
  root -b -q -l "$HERE/pp/ex_C15.C(\"$CSV\",\"$GFIT/\",${EBEAM}.0,5.0,\"_gated\",1.007825,15.0105993,\"15C(p,p') IC+PID gated\",\"ukf\")" \
       > "$LOG/ex_gated.log" 2>&1
  grep -E "protons -> Ex|good track" "$LOG/ex_gated.log" | tail -2
fi
echo "[$(date +%H:%M:%S)] GATED PASS DONE"
