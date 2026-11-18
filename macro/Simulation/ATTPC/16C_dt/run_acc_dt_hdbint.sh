#!/bin/bash
# 2026-11-18 repo copy of /mnt/f/a1975/sim_rebuild/dt_hdbint/run.sh (the copy that ran); outputs stay under that directory.
# 2026-11-18: (d,t) full-chain acceptance for the INTEGRAL-GATE production (gf_dt_hdbtau_int).
# Sims: the existing triton sims /mnt/f/a1975_C16_dt_sim/{gs,ex1..ex4}_s3001_sim.root (20000 ev each, NOT regenerated;
#       generated with the v2 gas volume, r 250 -- the data fit uses r280 geometry, as here).
# 1 RECO  like reco_d2_hdbtau: gain-100000 par, thr 40, HDBSCAN mcs 0 ms 3, float tau (as sim_rebuild/r280_hdbtau).
# 2 PLANE pid_width_dp.C on each sim reco -> same sqw (charge*tau_fit/4.5) as the data plane, to check the sim band vs the gate.
# 3 FIT   EXACTLY gf_dt_hdbtau_int.sh's arguments (r280, CATIMA, triton table, -2.85) twice:
#         gated  = Yassid's integral gate (frozen copy) + gateTauNomTB 4.5  -> fit_gate/
#         nogate = no PID gate                                                -> fit_nogate/
# 4 ACC   acceptance_fullchain_dt.C (chi2 < 4.5, vz 25-600 on truth z + 19.1, Ebeam 184.17) with an EMPTY gate on both fit
#         dirs: in fit_gate the gate already acted inside genfit, so its "fit" factor carries the gate.
# Resume = rerun. Each stage is read back before its marker is written.
source ~/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export ROOT_INCLUDE_PATH=$HOME/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
export GENFIT=/home/yassid/fair_install/GenFit; export LD_LIBRARY_PATH=$GENFIT/lib:${LD_LIBRARY_PATH:-}
export N=/mnt/f/a1975/sim_rebuild/dt_hdbint
export S=$HOME/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC/16C_dt
export D2=$HOME/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF
export GATE=/mnt/f/a1975/gf_dt_hdbtau_int/gate_used.json
mkdir -p $N/reco $N/fit_gate $N/fit_nogate $N/logs
exec 9>$N/.lock; flock -n 9 || { echo "another run.sh holds the lock"; exit 1; }
nent() { root -l -b -q -e "TFile f(\"$1\");auto*t=(TTree*)f.Get(\"cbmsim\");printf(\"NENT %lld\\n\",t?t->GetEntries():-1LL);" 2>/dev/null </dev/null | awk '/^NENT/{print $2}'; }
reco() { s=$1; SIM=/mnt/f/a1975_C16_dt_sim/${s}_s3001_sim.root; R=$N/reco/t_${s}_reco.root
  [ -e $R.done ] && { echo "have reco $s"; return 0; }
  cd $S; root -l -b -q "run_reco_C16dt.C(\"$SIM\",\"$R\",\"ATTPC.a1975_deuterium_dv1104_g100000.par\",40,0,3,0,20.0,\"mover\",0,\"hdbscan\",kTRUE,0.1)" > $N/logs/reco_$s.log 2>&1 </dev/null
  n0=$(nent $SIM); n=$(nent $R)
  grep -q "float-tau" $N/logs/reco_$s.log && [ -n "$n" ] && [ "$n" = "$n0" ] && { touch $R.done; echo "RECO $s ok $n $(date +%T)"; } || { echo "FAIL reco $s n=$n/$n0"; return 1; }; }
plane() { s=$1; R=$N/reco/t_${s}_reco.root; P=$N/reco/pid_width_t_${s}.root
  [ -e $P.done ] && return 0; cd $D2
  root -l -b -q "pid_width_dp.C(\"$R\",\"$P\")" > $N/logs/plane_$s.log 2>&1 </dev/null && grep -q tracks $N/logs/plane_$s.log && touch $P.done && echo "PLANE $s ok"; }
fit() { s=$1; mode=$2; R=$N/reco/; O=$N/fit_$mode/; F=${O}t_${s}_genfitter_t.root; L=$N/logs/fit_${mode}_$s.log
  [ -e $F.done ] && { echo "have fit $mode $s"; return 0; }
  if [ $mode = gate ]; then g="$GATE"; tn=4.5; else g=""; tn=0; fi
  cd $D2; root -l -b -q "fitGenfitter_a1975_deuterium.C(\"t_${s}\",-1,\"$R\",\"\",\"$O\",-2.85,2,5,\"$g\",4.0,10.0,170.0,kTRUE,kTRUE,1000010030,3.01550072,1,\"t\",\"_reco\",\"ATTPC_D300torr_r280_geomanager.root\",kTRUE,0,2,\"ATTPC.a1975_deuterium_dv1104.par\",kFALSE,kFALSE,kFALSE,\"triton_D2_300torr.txt\",kTRUE,kTRUE,kTRUE,kFALSE,0,kFALSE,0,kFALSE,5.0,25.0,0,95.0,8,kFALSE,kFALSE,\"\",kFALSE,$tn)" > $L 2>&1 </dev/null
  vr=$(root -b -l -q "pid/root_ok.C(\"$F\")" 2>/dev/null | grep -E '^(VALID|INVALID)' | head -1)
  ok=1; [[ "$vr" == VALID* ]] || ok=0; grep -q "Done\." $L || ok=0; grep -q "dE/dx from CATIMA" $L || ok=0
  grep -qi "segmentation violation" $L && ok=0; [ $mode = gate ] && { grep -q "TRACE-INTEGRAL" $L || ok=0; }
  [ $ok = 1 ] && { touch $F.done; echo "FIT $mode $s ok $vr $(date +%T)"; } || { echo "FAIL fit $mode $s ${vr:-unreadable}"; rm -f $F; return 1; }; }
export -f nent reco plane fit
printf '%s\n' gs ex1 ex2 ex3 ex4 | xargs -P 5 -I{} bash -c 'reco {}'
echo STAGE_RECO
printf '%s\n' gs ex1 ex2 ex3 ex4 | xargs -P 5 -I{} bash -c 'plane {}'
echo STAGE_PLANE
# 2026-11-18: the GATED fit is held (MODES=nogate by default): on 200 test events the sim integral band is 0.89-0.94 of the
# data band, so the data gate cuts 20-56 % of sim tritons. Rerun with MODES="gate nogate" (and the agreed gate handling).
MODES=${MODES:-nogate}
for s in gs ex1 ex2 ex3 ex4; do for m in $MODES; do echo "$s $m"; done; done | xargs -P 7 -L 1 bash -c 'fit "$@"' _
echo STAGE_FIT
cd $S
for m in $MODES; do
  root -l -b -q "acceptance_fullchain_dt.C(\"$N/reco/\",\"$N/fit_$m/\",\"\",\"$N/acc_dt_$m.root\")" > $N/acc_$m.txt 2>&1 </dev/null
  grep OVERALL $N/acc_$m.txt | sed "s/^/[$m] /"
done
echo ALLDONE
