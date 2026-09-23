#!/usr/bin/env bash
# Digitise + reconstruct the a1975 16C(d,p)17C simulation with the SAME chain the DATA uses.
#
# It calls run_reco_C16dt.C (the (d,t) sibling) unchanged -- that macro is fully parameterised and
# there is no (d,p)-specific reconstruction. Only the par file and the three chain settings differ.
#
# *** THE CHAIN SETTINGS ARE THE DATA'S, NOT THE (d,t) SIM'S ***
#   data production (D2_UKF/dv_dt.sh, reco_d2_dv1104):  psa multifit, thr 40, hdbscan mcs 0, ms 3
#   run_reco_C16dt.C DEFAULTS:                          thr 20, mcs 20, ms 8
# The (d,t) sim reco batch calls the macro with its defaults, so the (d,t) acceptance was measured
# on a chain that does NOT match the data it corrects -- the very mismatch run_reco_C16dt.C's own
# header warns about for (p,d). Not fixed here (it would invalidate the adopted (d,t) gain, which
# was calibrated at thr 20); flagged so it is not inherited silently. For (d,p) the data settings
# are passed explicitly.
#
# *** GAIN 75000, MEASURED FOR PROTONS 2026-09-23 -- DO NOT INHERIT THE (d,t) VALUE ***
# gain_check_dp.C on gs_s4001 (300 events, data-matched chain) vs run_0016, non-beam tracks,
# tracks >500 hits dropped, medians:
#       gain      sim hits/track   sim/data     clusters/track sim/data
#       35000          42.5          0.55              0.81     <- the (d,t) value: HALF the data
#       75000          80.0          1.04              1.32     <- ADOPTED
#      100000          94.0          1.22              1.39
#      150000         109.0          1.42              1.39     <- the (p,d) sim's value
# 2.1x the (d,t) gain, in the expected direction: a proton deposits less per pad than a triton and
# needs more gain to cross the same PSA threshold. It sits between the (d,t)'s 35000 and the (p,d)
# sim's 150000, both of which are now accounted for.
#
# WHAT THIS CALIBRATION DOES NOT ESTABLISH:
#   * Clusters/track is 1.32-1.39 at EVERY gain from 75k to 150k, i.e. flat -- so that mismatch is
#     NOT gain. The sim's clusters are less fragmented than the data's. Do not read a
#     cluster-count-dependent result off this sample without accounting for it.
#   * The medians match but the WIDTHS do not: at 75000 sim p25-p75 = 57-103 against data 32-169.
#     The sim is pure g.s. protons; the data is every non-beam track of every species. There is no
#     species selection on the data side without the PID gate. This is a calibration of chain
#     output, not a demonstration that the two samples are the same thing.
#   * It is ONE run (run_0016) and 300 simulated events.
#
#   ./run_reco_dp.sh [simFile] [gain] [nEvents]
set -eo pipefail
SIM=${1:-/mnt/f/a1975_C16_dp_sim/gs_s4001_sim.root}
GAIN=${2:-75000}
NEV=${3:-0}            # 0 = all
THR=40                 # data value
MCS=0                  # data value (0 = adaptive)
MS=3                   # data value
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../../../.." && pwd)
DT="$REPO/macro/Simulation/ATTPC/16C_dt"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"

PAR="ATTPC.a1975_deuterium_dv1104_g${GAIN}.par"
[ -s "$REPO/parameters/$PAR" ] || { echo "ERROR: no $PAR -- run 16C_dt/make_gain_pars.sh $GAIN"; exit 1; }
[ -s "$SIM" ] || { echo "ERROR: no simulation at $SIM"; exit 1; }

OUT="${SIM%_sim.root}_reco_g${GAIN}.root"
LOG="${SIM%_sim.root}_reco_g${GAIN}.log"
MARK="${OUT}.done"

if [ -f "$MARK" ]; then echo "[$(date +%H:%M:%S)] already done: $OUT"; exit 0; fi

echo "[$(date +%H:%M:%S)] reco  $(basename "$SIM")  gain=$GAIN  thr=$THR mcs=$MCS ms=$MS  -> $(basename "$OUT")"
cd "$DT"
# ROOT's exit code is not a verdict -- it segfaults in EndOfProcessCleanups after writing a good
# file (that is how a complete 20000-event simulation got reported as FAILED). The readback below
# is the verdict. See run_gs_dp.sh for the same trap.
root -b -l -q "run_reco_C16dt.C(\"$SIM\",\"$OUT\",\"$PAR\",$THR,$MCS,$MS,$NEV)" > "$LOG" 2>&1 || true

[ -s "$OUT" ] || { echo "[$(date +%H:%M:%S)] FAILED: no output; see $LOG"; exit 1; }
# *** THE READBACK NEEDS THE SAME `|| true` AS THE MAIN CALL ***
# I guarded the reco call and not its verification, and the verification root ALSO segfaults at
# exit. Under `set -eo pipefail` the failed assignment killed the script, so a complete, valid
# 8.18 GB reconstruction of 20000 events was reported FAILED (exit 17) -- the second time the same
# trap bit in one session. Anything that shells out to root under `set -e` needs this.
RB=$(root -b -l -q -e "TFile*F=TFile::Open(\"$OUT\");TTree*t=(TTree*)F->Get(\"cbmsim\");printf(\"NEVT %lld REC %d\n\",t?t->GetEntries():-1,F->TestBit(TFile::kRecovered));" 2>/dev/null || true)
N=$(printf '%s' "$RB" | grep -oP 'NEVT \K-?\d+' || true)
R=$(printf '%s' "$RB" | grep -oP 'REC \K\d+' || true)
if [ "${N:--1}" -lt 1 ] || [ "${R:-1}" -ne 0 ]; then
  echo "[$(date +%H:%M:%S)] FAILED readback (entries=$N recovered=$R); deleting"; rm -f "$OUT"; exit 1
fi
touch "$MARK"
echo "[$(date +%H:%M:%S)] DONE  $N events  $(du -h "$OUT" | cut -f1)  -> $OUT"
