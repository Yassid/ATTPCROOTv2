#!/bin/bash
# 2026-11-16: the r280 (d,p) sims (sim_rebuild/r280/p_gs_sim.root, p_ex39_sim.root -- NOT regenerated) re-reconstructed
# like the HDBSCAN + float-tau D2 data (praType "hdbscan" mcs 0 ms 3, floatTau; TC dropped 2026-11-16: it fragments backward spirals), gain 100000 par, thr 40; then EXACTLY
# r280/run_sim.sh's fit (-2.85, no gate, FIND 80->5, r280 geometry) and sim_dump_dp.C. Each stage read back; resumable.
source /mnt/f/a1975/sim_rebuild/r280/env.sh; S=$HOME/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC
N=/mnt/f/a1975/sim_rebuild/r280_hdbtau
nent() { root -l -b -q -e "TFile f(\"$1\");auto*t=(TTree*)f.Get(\"cbmsim\");printf(\"NENT %lld\\n\",t?t->GetEntries():-1LL);" 2>/dev/null </dev/null | awk '/^NENT/{print $2}'; }
one() {
  J=$1; SIM=/mnt/f/a1975/sim_rebuild/r280/${J}_sim.root; RECO=$N/${J}_reco.root; F=$N/${J}_genfitter_p.root
  n0=$(nent $SIM)
  if [ ! -e $RECO.done ]; then cd $S/16C_dt; root -l -b -q "run_reco_C16dt.C(\"$SIM\",\"$RECO\",\"ATTPC.a1975_deuterium_dv1104_g100000.par\",40,0,3,0,20.0,\"mover\",0,\"hdbscan\",kTRUE,0.1)" > $N/logs/reco_$J.log 2>&1 </dev/null
    n=$(nent $RECO); ! grep -q "AtTrackFinderTC" $N/logs/reco_$J.log && grep -q "float-tau" $N/logs/reco_$J.log && [ "$n" = "$n0" ] && { touch $RECO.done; echo "RECO $J ok $n $(date +%T)"; } || { echo "FAIL reco $J n=$n/$n0"; return 1; }; fi
  if [ ! -e $F.done ]; then cd $D2; root -l -b -q "fitGenfitter_a1975_deuterium.C(\"$J\",-1,\"$N/\",\"\",\"$N/\",-2.85,2,5,\"\",4.0,10.0,170.0,kTRUE,kTRUE,2212,1.00782503207,1,\"p\",\"_reco\",\"ATTPC_D300torr_r280_geomanager.root\",kTRUE,0,2,\"ATTPC.a1975_deuterium_dv1104.par\",kFALSE,kFALSE,kFALSE,\"proton_D2_300torr.txt\",kTRUE,kTRUE,kTRUE,kFALSE,0,kTRUE,0.1,kFALSE,5.0,5.0,1,80.0,8)" > $N/logs/fit_$J.log 2>&1 </dev/null
    vr=$(root -b -l -q "pid/root_ok.C(\"$F\")" 2>/dev/null | grep -E '^(VALID|INVALID)' | head -1); [[ "$vr" == VALID* ]] && { touch $F.done; echo "FIT $J ok $vr $(date +%T)"; } || { echo "FAIL fit $J ${vr:-unreadable}"; return 1; }; fi
  cd $S/16C_dp; root -l -b -q "sim_dump_dp.C(\"$F\",\"$SIM\",\"$N/${J}_dump.root\",true)" </dev/null > $N/logs/dump_$J.log 2>&1; tail -1 $N/logs/dump_$J.log; echo "DUMP $J $(date +%T)"
}
export -f one nent; export N S
printf '%s\n' p_gs p_ex39 | xargs -P 2 -I{} bash -c 'one {}'
echo SIMHDBTAU_ALLDONE
