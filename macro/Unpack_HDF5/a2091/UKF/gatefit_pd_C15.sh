#!/bin/bash
# Gated production for a2091 15C(p,d)14C: IC(15C) + DEUTERON PID gate -> UKF(deuteron) -> Ex.
#
# Same proton-target runs and the same reco as the (p,p') pass -- only the ejectile differs,
# so nothing is re-unpacked. Outputs go to their own directories so the (p,p') results are
# never touched.
#
# This replaces fitpipe_pd_C15.sh, which carried three stale assumptions from the port:
#   - it read a2091_C15_reco_slim/, a slim cache that was NEVER built for a2091 (every run
#     would have hit "SKIP (missing)" and produced nothing, silently)
#   - its IC window was 500-900, not the 979.5-1278.8 window measured for the 15C beam
#   - genfit_pd_C15.sh defaulted to ATTPC_H600torr geometry; a2091 ran at 300 torr
#
# Masses use the ATOMIC values, consistent with the (p,p') pass using 1.007825 for the proton
# and 15.0105993 for 15C:  deuteron 2.014102, 14C 14.003242. (ex_C15.C's docstring suggests
# 2.013553, which is the deuteron NUCLEAR mass -- mixing conventions would shift Ex.)
set -u
REPO=/home/yassid/fair_install/ATTPCROOTv2
HERE=$REPO/macro/Unpack_HDF5/a2091/UKF
RECO=/home/yassid/a2091_C15_reco            # symlink -> Seagate, shared with (p,p')
IC=/home/yassid/a2091_C15_ic
GIN=/home/yassid/a2091_C15_fit/in_pd        # gated (p,d) recos
GFIT=/home/yassid/a2091_C15_fit/gated_pd    # gated (p,d) fits
LOG=/home/yassid/a2091_C15_fit/logs/gated_pd
DGATE=$HERE/pid/deuteron_15C.json
NPAR="${2:-6}"
EBEAM="${EBEAM:-170}"                       # same beam as (p,p'); (p,d) is an independent check
DENS="${DENS:-3.308e-5}"
THMIN="${THMIN:-90}"                        # Spyral polar is the 180-theta mirror, so forward
                                            # deuterons sit at HIGH polar and pass this
MEJ=2.014102                                # deuteron (atomic)
MRES=14.003242                              # 14C   (atomic)
MINRECO=5000000
mkdir -p "$GIN" "$GFIT" "$LOG"

trap 'pkill -P $$ 2>/dev/null' EXIT INT TERM

set +u
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
source "$REPO/build/config.sh" >/dev/null 2>&1
export ROOT_INCLUDE_PATH="$REPO/build/include:/home/yassid/fair_install/FairRootInstall/include"
set -u

[ -f "$DGATE" ] || { echo "FATAL: no deuteron gate at $DGATE"; exit 1; }

if [ -n "${1:-}" ]; then RUNS="$1"; else
  RUNS=$(for f in "$RECO"/*_reco.root; do
           [ -f "$f" ] || continue
           [ "$(stat -c%s "$f")" -ge "$MINRECO" ] || continue
           r=$(basename "$f" _reco.root)
           [ -f "$IC/${r}_ic.root" ] && echo "$r"
         done | sort -u | tr '\n' ' ')
fi
echo "[$(date +%H:%M:%S)] (p,d) gated pass over $(echo $RUNS | wc -w) runs, NPAR=$NPAR, Ebeam=$EBEAM, thMin=$THMIN"

one(){ local r="$1"
  if [ ! -f "$GIN/${r}_reco.root" ]; then
    root -b -q -l "$HERE/pipeline/gate_events_C15.C(\"$r\",\"$RECO/\",\"$GIN/\",\"$RECO/${r}_reco.root\",979.5,1278.8,1050,1250,2.85,\"$DGATE\",$THMIN)" \
         > "$LOG/${r}_gate.log" 2>&1
  fi
  [ -f "$GIN/${r}_reco.root" ] || { echo "[$(date +%H:%M:%S)] $r GATE FAILED"; return; }
  if [ ! -f "$GFIT/${r}_ukf.root" ]; then
    root -b -q -l "$HERE/pipeline/fitUKF_C15.C(\"$r\",-1,\"deuteron\",-1,2.85,$DENS,\"\",\"$GIN/\",0.5,0.1,1,10,\"$GFIT/\")" \
         > "$LOG/${r}_ukf.log" 2>&1
  fi
  echo "[$(date +%H:%M:%S)] $r $(grep -oE 'gated-proton events, [0-9]+ tracks' "$LOG/${r}_gate.log" 2>/dev/null | head -1) ukf=$([ -f "$GFIT/${r}_ukf.root" ] && echo ok || echo FAIL)"
}
export -f one; export HERE RECO IC GIN GFIT LOG DENS DGATE THMIN
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}

echo "[$(date +%H:%M:%S)] (p,d) gated: $(ls "$GIN"/*_reco.root 2>/dev/null | wc -l) recos, $(ls "$GFIT"/*_ukf.root 2>/dev/null | wc -l) ukf fits"

CSV=$(ls "$GFIT"/*_ukf.root 2>/dev/null | sed 's#.*/##;s/_ukf.root//' | sort -u | paste -sd,)
if [ -n "$CSV" ]; then
  echo "[$(date +%H:%M:%S)] Ex spectrum, (p,d), Ebeam=$EBEAM"
  root -b -q -l "$HERE/pp/ex_C15.C(\"$CSV\",\"$GFIT/\",${EBEAM}.0,5.0,\"_pd_ukf\",$MEJ,$MRES,\"15C(p,d)14C\",\"ukf\")" \
       > "$LOG/ex_pd.log" 2>&1
  grep -E "protons -> Ex|good track|channel:" "$LOG/ex_pd.log" | tail -3
fi
echo "[$(date +%H:%M:%S)] (p,d) PASS DONE"
