#!/bin/bash
# 16C(d,d) elastic simulation for the scale-test acceptance, r280 geometry (2026-11-10).
# gen C16_dd_el_sim_r280.C (theta_cm 8-40, Ebeam 184.25 as p_gs r280) -> reco at gain 100000 (run_sim.sh settings) ->
# deuteron genfit with EXACTLY the arguments of run_prod_dd.sh (data), no gate. Each stage read back; resumable via .done.
# Usage: run_dd_el.sh OUTDIR NEV "SEEDS" NPAR
source /mnt/f/a1975/sim_rebuild/r280/env.sh
S=$HOME/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC
OUT=${1:-$O/dd_el}; NEV=${2:-5000}; SEEDS=${3:-"6001 6002 6003 6004"}; NPAR=${4:-2}
mkdir -p $OUT/logs
nent() { root -l -b -q -e "TFile f(\"$1\");auto*t=(TTree*)f.Get(\"cbmsim\");printf(\"NENT %lld\\n\",t?t->GetEntries():-1LL);" 2>/dev/null </dev/null | awk '/^NENT/{print $2}'; }
one() {
  s=$1; SIM=$OUT/el_s${s}_sim.root; RECO=$OUT/el_s${s}_reco.root; F=$OUT/el_s${s}_genfitter_d.root
  if [ ! -e $SIM.done ]; then cd $S/16C_dt
    root -b -l -q "C16_dd_el_sim_r280.C($NEV,8.0,40.0,\"TGeant4\",-28.5,\"$SIM\",0.0,$s,184.25)" > $OUT/logs/sim_$s.log 2>&1 </dev/null
    n=$(nent $SIM); grep -q "RNG seed requested: $s" $OUT/logs/sim_$s.log && [ "$n" = "$NEV" ] && { touch $SIM.done; echo "SIM $s ok $n $(date +%T)"; } || { echo "FAIL sim $s n=$n"; return; }
  fi
  if [ ! -e $RECO.done ]; then cd $S/16C_dt
    root -l -b -q "run_reco_C16dt.C(\"$SIM\",\"$RECO\",\"ATTPC.a1975_deuterium_dv1104_g100000.par\",40,0,3)" > $OUT/logs/reco_$s.log 2>&1 </dev/null
    n=$(nent $RECO); grep -q "par = ATTPC.a1975_deuterium_dv1104_g100000.par" $OUT/logs/reco_$s.log && [ "$n" = "$NEV" ] && { touch $RECO.done; echo "RECO $s ok $n $(date +%T)"; } || { echo "FAIL reco $s n=$n"; return; }
  fi
  if [ ! -e $F.done ]; then cd $D2
    root -l -b -q "fitGenfitter_a1975_deuterium.C(\"el_s${s}\",-1,\"$OUT/\",\"\",\"$OUT/\",-2.85,2,5,\"\",4.0,10.0,170.0,kTRUE,kTRUE,1000010020,2.01355321,1,\"d\",\"_reco\",\"ATTPC_D300torr_r280_geomanager.root\",kTRUE,0,2,\"ATTPC.a1975_deuterium_dv1104.par\",kFALSE,kFALSE,kFALSE,\"deuteron_D2_300torr.txt\",kTRUE,kTRUE,kTRUE,kFALSE)" > $OUT/logs/fit_$s.log 2>&1 </dev/null
    vr=$(root -b -l -q "pid/root_ok.C(\"$F\")" 2>/dev/null | grep -E '^(VALID|INVALID)' | head -1)
    [[ "$vr" == VALID* ]] && grep -q "dE/dx from CATIMA" $OUT/logs/fit_$s.log && { touch $F.done; echo "FIT $s ok $vr $(date +%T)"; } || { echo "FAIL fit $s ${vr:-unreadable}"; return; }
  fi
  echo "DONE $s $(date +%T)"
}
export -f one nent; export OUT NEV S
printf '%s\n' $SEEDS | xargs -P $NPAR -I{} bash -c 'one {}'
echo DD_EL_ALLDONE $(ls $OUT/*_genfitter_d.root.done 2>/dev/null | wc -l)
