#!/bin/bash
source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh 2>/dev/null
cd /home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
for sig in 0 1.0 1.5 2.0 2.5; do
  pfx="s$(echo $sig | tr -d '.')_"
  root -l -b -q "unpackReco_multifit.C(\"run_0305\", 200, false, \"/tmp/claude-1000/-home-yassid/b5b91330-d765-4197-9f1c-6da8c3d5a215/scratchpad/$pfx\", \"/home/yassid/spyral_d2/h5/\", false, true, \"multifit\", $sig)" > /dev/null 2>&1
  echo "done sigma=$sig"
done
echo "SCAN_DONE"
