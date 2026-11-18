#!/bin/bash
# 2026-11-18 repo copy of /mnt/f/a1975/sim_rebuild/dt_hdbint/run_gate.sh (the copy that ran); outputs stay under that directory.
# 2026-11-18: GATED sim fit, sim ALIGNED to data (Yassid: "align the sim to data and apply the gate"; "one global factor if
# the slices agree"). Sim/data triton-band median ratio on the integral plane, Brho 0.30-1.40 (9 slices, full sim planes):
# .914 .927 .929 .915 .904 .866 .860 .873 .906, weighted mean 0.886 (spread +-4.8 %) -> gateTauNomTB = 4.5 * 0.886^2 = 3.534,
# i.e. the gate sees sim charge x tau_fit / 3.534 (sqrt(dE/dx) up by 1/0.886), Brho untouched.
# Same fitter arguments as run.sh / gf_dt_hdbtau_int.sh otherwise. Then acceptance_fullchain_dt.C with an empty gate.
source ~/fair_install/ATTPCROOTv2-OpenKF/build/config.sh >/dev/null 2>&1
export ROOT_INCLUDE_PATH=$HOME/fair_install/ATTPCROOTv2-OpenKF/build/include:${ROOT_INCLUDE_PATH:-}
export GENFIT=/home/yassid/fair_install/GenFit; export LD_LIBRARY_PATH=$GENFIT/lib:${LD_LIBRARY_PATH:-}
export N=/mnt/f/a1975/sim_rebuild/dt_hdbint S=$HOME/fair_install/ATTPCROOTv2-OpenKF/macro/Simulation/ATTPC/16C_dt
export D2=$HOME/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF GATE=/mnt/f/a1975/gf_dt_hdbtau_int/gate_used.json TN=3.534
O=$N/fit_gate_aligned; mkdir -p $O; export O
exec 9>$N/.lock_gate; flock -n 9 || { echo "another run_gate.sh holds the lock"; exit 1; }
fitg() { s=$1; F=$O/t_${s}_genfitter_t.root; L=$N/logs/fit_gate_aligned_$s.log
  [ -e $F.done ] && { echo "have $s"; return 0; }
  cd $D2; root -l -b -q "fitGenfitter_a1975_deuterium.C(\"t_${s}\",-1,\"$N/reco/\",\"\",\"$O/\",-2.85,2,5,\"$GATE\",4.0,10.0,170.0,kTRUE,kTRUE,1000010030,3.01550072,1,\"t\",\"_reco\",\"ATTPC_D300torr_r280_geomanager.root\",kTRUE,0,2,\"ATTPC.a1975_deuterium_dv1104.par\",kFALSE,kFALSE,kFALSE,\"triton_D2_300torr.txt\",kTRUE,kTRUE,kTRUE,kFALSE,0,kFALSE,0,kFALSE,5.0,25.0,0,95.0,8,kFALSE,kFALSE,\"\",kFALSE,$TN)" > $L 2>&1 </dev/null
  vr=$(root -b -l -q "pid/root_ok.C(\"$F\")" 2>/dev/null | grep -E '^(VALID|INVALID)' | head -1)
  ok=1; [[ "$vr" == VALID* ]] || ok=0; for k in "Done\." "dE/dx from CATIMA" "TRACE-INTEGRAL"; do grep -q "$k" $L || ok=0; done
  grep -qi "segmentation violation" $L && ok=0
  [ $ok = 1 ] && { touch $F.done; echo "FIT $s ok $vr $(date +%T)"; } || { echo "FAIL fit $s ${vr:-unreadable}"; rm -f $F; return 1; }; }
export -f fitg
printf '%s\n' gs ex1 ex2 ex3 ex4 | xargs -P 5 -I{} bash -c 'fitg {}'
echo STAGE_FIT
cd $S; root -l -b -q "acceptance_fullchain_dt.C(\"$N/reco/\",\"$O/\",\"\",\"$N/acc_dt_gate_aligned.root\")" > $N/acc_gate_aligned.txt 2>&1 </dev/null
grep OVERALL $N/acc_gate_aligned.txt; echo ALLDONE
