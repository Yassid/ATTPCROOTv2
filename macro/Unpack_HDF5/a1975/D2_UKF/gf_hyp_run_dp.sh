#!/bin/bash
# gf_hyp_run_dp.sh <ioDir> <fileName> <nEvents> <hyp d|t> <outDir> [measSigma, default 4.0]
# Genfit fit of a simulated reco file under ONE hypothesis, settings of the (d,t) production except CATIMA eloss over
# the full range for both species and no PID gate (2026-11-01 two-hypothesis test). Example, all 12 fits:
#   for f in el_s5001..4 (/mnt/f/a1975_C16_dd_sim/), gs_s3001 ex2_s3001 (/mnt/f/a1975_C16_dt_sim/); h in d t
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export GENFIT=/home/yassid/fair_install/GenFit
export LD_LIBRARY_PATH=$GENFIT/lib:${LD_LIBRARY_PATH:-}
export ROOT_INCLUDE_PATH=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
IO=$1; FN=$2; NE=$3; H=$4; OUT=$5; SIG=${6:-4.0}
if [ "$H" = t ]; then PDG=1000010030; M=3.01550072; else PDG=1000010020; M=2.01355321; fi
root -l -b -q "fitGenfitter_a1975_deuterium.C(\"$FN\",$NE,\"$IO\",\"\",\"$OUT\",\
-2.85,2,5,\"\",$SIG,10.0,170.0,kTRUE,kTRUE,$PDG,$M,1,\"$H\",\"_reco\",\
\"ATTPC_D300torr_v2_geomanager.root\",kTRUE,0,2,\"ATTPC.a1975_deuterium_dv1104_g35000.par\",kFALSE,kFALSE,kFALSE,\"\",kTRUE,kTRUE,kTRUE,kTRUE)" < /dev/null
