#!/usr/bin/env bash
# Generate the a1975 16C(d,p)17C GROUND-STATE simulation.
#
# WHY IT EXISTS: the data's BOUND region shows a low-Ex peak that walks +0.30 MeV across the
# chamber and a resolution collapse at v_z 650-750 mm. Neither can be attributed to geometry or to
# reconstruction without truth carried through the same chain. See C16_dp_sim.C's header.
#
# THE CM RANGE IS DELIBERATELY NARROW -- 2 to 50 deg, chosen by Yassid 2026-09-23.
# Measured on a 600-event full-range smoke test (truth, kine2b, analysis convention):
#       theta_lab   theta_cm      proton KE
#         0- 60     80-171        37.9 MeV
#        60- 89     45- 80        15.2
#        89-110     30- 45         6.9     <-- the analysis window starts here
#       110-140     13- 28         3.4
#       140-161      9             1.8
# The adopted selection is theta_lab 89-161, so the whole analysis lives at theta_cm <= 50 and the
# full 2-178 range spends ~85 % of the CPU on forward-lab protons that are never used (only 14 % of
# generated protons landed in 89-161).
#
# *** WHAT THIS RANGE COSTS, and it must not be forgotten ***
#   1. The TOTAL ANGULAR acceptance cannot be measured from this sample -- only the acceptance and
#      resolution VS VERTEX Z inside the analysis window. That is the question it was built for.
#   2. Protons MIGRATING INTO the window from theta_cm > 50 via reconstruction error are not
#      generated, so any in/out migration is one-sided here.
#   The (d,t) sibling deliberately generates 2-178 for exactly these reasons; if a full acceptance
#   is ever wanted, regenerate at 2-178 rather than rescaling this one.
#
# CONVENTION: AtTPC2Body ranges over the RESIDUAL's cm angle, which is the same convention the
# analysis uses (theta_cm = pi - theta3_cm). Verified on the smoke test above -- generating 2-178
# returned analysis theta_cm 2.5-171.4.
#
#   ./run_gs_dp.sh [nEvents] [seed] [outDir]
set -eo pipefail
NEV=${1:-20000}
SEED=${2:-4001}
OUT=${3:-/mnt/f/a1975_C16_dp_sim}
EBEAM=184.25
CMLO=2.0
CMHI=50.0
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../../../.." && pwd)
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
mkdir -p "$OUT"
cd "$HERE"

J="gs_s${SEED}"
SIM="$OUT/${J}_sim.root"
LOG="$OUT/${J}_sim.log"

# A completion marker must be EARNED by reading the product back, never by the process exiting --
# the (d,p) production driver was marking truncated files [ok] for exactly this reason.
if [ -f "$OUT/$J.marker" ]; then
  echo "[$(date +%H:%M:%S)] $J already COMPLETED"; exit 0
fi

echo "[$(date +%H:%M:%S)] 16C(d,p)17C g.s.  N=$NEV  seed=$SEED  Ebeam=$EBEAM  CM=$CMLO-$CMHI"
echo "                    -> $SIM"
# *** ROOT'S EXIT CODE IS NOT A VERDICT HERE -- `|| true` IS DELIBERATE ***
# This macro finishes its 20000 events, writes the file, prints "Macro finished successfully",
# and THEN segfaults in TROOT::EndOfProcessCleanups, returning 11. Under `set -e` that killed the
# script before the readback below and reported a complete, valid 565 MB simulation as FAILED.
# Same shape as the `root exits 8 under set -e` trap in reference_campaign_pipeline_traps.
# The readback is the verdict -- that is the whole point of having one.
root -b -l -q "C16_dp_sim.C($NEV,$CMLO,$CMHI,\"TGeant4\",-28.5,\"$SIM\",0.0,$SEED,$EBEAM)" > "$LOG" 2>&1 || true
rc=$?
[ $rc -ne 0 ] && echo "[$(date +%H:%M:%S)] note: root exited $rc -- judging on the readback, not on this"

if [ ! -s "$SIM" ]; then
  echo "[$(date +%H:%M:%S)] FAILED: no output written; see $LOG"; exit 1
fi
# read the product back
# The readback root segfaults at exit too -- without `|| true` here, `set -eo pipefail` kills the
# script on a SUCCESSFUL verification. Same trap as the main call above.
N=$(root -b -l -q -e "TFile*F=TFile::Open(\"$SIM\");TTree*t=(TTree*)F->Get(\"cbmsim\");printf(\"NEVT %lld\n\",t?t->GetEntries():-1);" 2>/dev/null | grep -oP 'NEVT \K-?\d+' || true)
if [ "${N:--1}" -lt 1 ]; then
  echo "[$(date +%H:%M:%S)] FAILED readback (entries=$N); deleting"; rm -f "$SIM"; exit 1
fi
touch "$OUT/$J.marker"
echo "[$(date +%H:%M:%S)] DONE  $J  $N events  $(du -h "$SIM" | cut -f1)"
