#!/bin/bash
# ONE command to produce a simulated sample that looks like the recorded data.
#
#     ./run_simulation.sh            6000 events (3000 reactions), ~15 min
#     ./run_simulation.sh 30000      the high-statistics sample, ~50 min
#     ./run_simulation.sh 2000 out2  write somewhere other than ./sim_out
#
#     SETTING=b057 ./run_simulation.sh 6000 out_b057
#         the magnet-on setting: 150 torr, B = 0.5691 T.  Tracks are HELICES there, so the
#         straight-line finder in step 2 does not apply -- see the guide.
#
# Every physics parameter in here is CALIBRATED -- beam energy and tilt, gas, drift,
# PSA threshold, pile-up rate.  You should not need to touch any of them, and the guide
# says what each one is and how it was pinned down before you do.
#
# Produces, in the output directory:
#     sim.root      Geant4 truth: what actually happened
#     merged.root   beam + two arms per event, no pile-up
#     real.root     the same with pile-up added -- USE THIS ONE
#     truth.csv     per-reaction E_cm, theta_cm, directions and vertex
#
# Then measure the tracking efficiency with, from the analysis package:
#     ./sim_efficiency.sh <outdir>/real.root <outdir>/truth.csv
# Source the ATTPCROOT environment BEFORE "set -u".  setup_fr19port.sh expands variables
# that may be unset, and with -u in force that kills the shell -- silently, because its
# stderr is redirected here.  The symptom is a script that prints nothing and exits 0.
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
source /home/yassid/fair_install/ATTPCROOTv2_fr19port/setup_fr19port.sh >/dev/null 2>&1
set -u

N=${1:-6000}
D=${2:-sim_out}
PKG=${DEC2014_PKG:-/home/yassid/dec2014_package}
PILEUP=${PILEUP:-0.10}        # measured from the data; see the guide
SETTING=${SETTING:-b0}        # b0 = 299.5 torr B=0;  b057 = 150 torr B=0.5691 T
case "$SETTING" in
  b0)   GEN=He4He4_sim_el_11p4.C; DIG=rundigi_sim_11p4.C; GENDIR=data11p4 ;;
  b057) GEN=He4He4_sim_b057.C;    DIG=rundigi_sim_b057.C; GENDIR=data_b057 ;;
  *) echo "unknown SETTING=$SETTING (expected b0 or b057)"; exit 1 ;;
esac
mkdir -p "$GENDIR"
mkdir -p "$D"

say() { echo "[$(date +%T)] $*"; }

say "1/5  Geant4 [$SETTING]: $N events (half are beam traverses, half are reactions)"
root -l -b -q "$GEN($N)" > "$D/sim.log" 2>&1 \
  || { echo "  FAILED -- see $D/sim.log"; exit 1; }
cp "$GENDIR/attpcsim_in.root" "$D/sim.root"

say "2/5  digitising (calibrated: PSA threshold 1000, bulk branches off)"
root -l -b -q "$DIG(\"./$D/sim.root\",\"./$D/digi.root\")" > "$D/digi.log" 2>&1 \
  || { echo "  FAILED -- see $D/digi.log"; exit 1; }

say "3/5  merging each beam traverse with its reaction"
root -l -b -q "merge_beam_reaction.C(\"$D/digi.root\",\"$D/merged.root\")" > "$D/merge.log" 2>&1
grep -E "merged events|hits/event" "$D/merge.log" | sed 's/^/       /'

say "4/5  adding beam pile-up at the measured rate ($PILEUP)"
root -l -b -q "add_pileup.C(\"$D/merged.root\",\"$D/digi.root\",\"$D/real.root\",$PILEUP)" 2>&1 \
  | grep -E "hits/event|events," | sed 's/^/       /'

say "5/5  extracting MC truth"
root -l -b -q "$PKG/lib/mc_truth2.C(\"$D/sim.root\",\"$D/truth.csv\")" 2>&1 \
  | grep -E "reactions" | sed 's/^/       /'

say "DONE.  Analyse $D/real.root against $D/truth.csv"
du -sh "$D"/*.root "$D"/truth.csv 2>/dev/null | sed 's/^/       /'
