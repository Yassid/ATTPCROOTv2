#!/usr/bin/env bash
# 2026-11-16 11:55: defaults point to the HDBSCAN + float-tau products (TC dropped for (d,p)); directory name kept.
# 2026-11-16: inputs for the (d,p) bound AD under TC + float-tau conditions, same recipe as ../trig (dump_keys.C):
#   gf_data.txt = mkexp_pp.C(cache, IC 900-1300) slim -> ke th vz c2   (cache = caches/dp_kin_hdbtau.root)
#   gf_sim.txt  = r280_hdbtau/p_gs_dump.root -> ke th vz c2 tke tth tvz
#   keys.txt    = copied from ../trig (same truth: sims are NOT regenerated, only re-reconstructed)
set +u; source /home/yassid/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1; set -u
D=${D:-/mnt/f/a1975/sim_rebuild/trig_tctau}; PP=/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/UKF/pp
CACHE=${CACHE:-/mnt/f/a1975/caches/dp_kin_hdbtau.root}; SDUMP=${SDUMP:-/mnt/f/a1975/sim_rebuild/r280_hdbtau/p_gs_dump.root}
if [ -s "$CACHE" ]; then
  root -b -l -q "$PP/mkexp_pp.C(\"$CACHE\",\"$D/exp_dp_tctau.root\",1e9,900,1300)" > $D/mkexp.log 2>&1
  root -b -l -q -e "TFile f(\"$D/exp_dp_tctau.root\"); auto*t=(TTree*)f.Get(\"pk\"); float ke,th,vz,c2; t->SetBranchAddress(\"ke\",&ke); t->SetBranchAddress(\"theta\",&th); t->SetBranchAddress(\"vertexz\",&vz); t->SetBranchAddress(\"chi2ndf\",&c2); FILE*o=fopen(\"$D/gf_data.txt\",\"w\"); for(long i=0;i<t->GetEntries();++i){t->GetEntry(i); fprintf(o,\"%.3f %.2f %.1f %.3f\n\",ke,th,vz,c2);} fclose(o);"
  echo "gf_data $(wc -l < $D/gf_data.txt)"
fi
if [ -s "$SDUMP" ]; then
  root -b -l -q -e "TFile f(\"$SDUMP\"); auto*t=(TTree*)f.Get(\"r\"); double v[7]; const char*N[7]={\"ke\",\"th\",\"vz\",\"c2\",\"tke\",\"tth\",\"tvz\"}; for(int i=0;i<7;++i) t->SetBranchAddress(N[i],&v[i]); FILE*o=fopen(\"$D/gf_sim.txt\",\"w\"); for(long i=0;i<t->GetEntries();++i){t->GetEntry(i); fprintf(o,\"%.4f %.3f %.2f %.4f %.4f %.4f %.4f\n\",v[0],v[1],v[2],v[3],v[4],v[5],v[6]);} fclose(o);"
  echo "gf_sim $(wc -l < $D/gf_sim.txt)"
fi
