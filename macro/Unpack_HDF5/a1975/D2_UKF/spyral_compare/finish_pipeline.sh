#!/usr/bin/env bash
# Wait for the loose PRA reco to finish, then: fit it, build Ex ntuples, compute
# spiral flags, and produce the default-vs-loose-vs-adaptive Ex comparison.
set -e
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
IODIR=/mnt/f/a1975/reco_d2/
D=$REPO/macro/Unpack_HDF5/a1975/D2_UKF
SC=$D/spyral_compare
RECO_PID=68908

echo "[wait] loose reco (PID $RECO_PID) ... $(date '+%H:%M:%S')"
while kill -0 $RECO_PID 2>/dev/null; do sleep 30; done
touch ${IODIR}run_0016loose_reco.root.done
echo "[ok] reco done $(date '+%H:%M:%S')"

cd $REPO/build && source config.sh >/dev/null 2>&1
export LD_LIBRARY_PATH="/home/yassid/fair_install/GenFit/lib:$LD_LIBRARY_PATH"
cd $D

echo "[fit] loose genfit $(date '+%H:%M:%S')"
root -l -b -q "fitGenfitter_a1975_deuterium.C(\"run_0016loose\", -1, \"$IODIR\", \"\", \"$IODIR\", -2.85)" \
   > $SC/fit_loose.log 2>&1
echo "[fit] done $(date '+%H:%M:%S')"

echo "[ex] build loose Ex ntuple"
root -l -b -q 'spyral_compare/scripts/build_ex_ntuple.C("run_0016loose","","spyral_compare/cache_loose.root")' \
   >> $SC/fit_loose.log 2>&1

echo "[flag] spiral flags"
root -l -b -q 'spyral_compare/scripts/spiral_flag.C("run_0016")' >> $SC/fit_loose.log 2>&1

echo "[cmp] Ex comparison"
~/Spyral/venv/bin/python spyral_compare/scripts/compare_ex.py > $SC/compare_ex.out 2>&1

echo "===== PIPELINE DONE $(date '+%H:%M:%S') ====="
cat $SC/compare_ex.out
