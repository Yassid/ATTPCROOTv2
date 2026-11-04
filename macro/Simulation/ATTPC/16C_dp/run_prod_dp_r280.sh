#!/bin/bash
# All find805 data runs refitted with gas to r 280 mm instead of the v2 steel wall at r 250. Same arguments as prod_dp_findscan.sh find805; only the geometry and the output dir differ. New dir, /mnt/f/a1975/gf_dp_find805 untouched. 2026-11-04.
source /mnt/f/a1975/sim_rebuild/r280/env.sh; cd $D2
G=/mnt/f/a1975/gf_dp_find805_r280
one() { r=$1; f=$G/${r}_multifit_genfitter_p.root; if [ -e $f.done ]; then echo "have $r"; return; fi
  root -l -b -q "fitGenfitter_a1975_deuterium.C(\"${r}_multifit\",-1,\"/mnt/f/a1975/reco_d2_dv1104/\",\"\",\"$G/\",-2.85,2,5,\"pid/proton_d2_dv1104.json\",4.0,10.0,170.0,kTRUE,kTRUE,2212,1.00782503207,1,\"p\",\"_reco\",\"ATTPC_D300torr_r280_geomanager.root\",kTRUE,0,2,\"ATTPC.a1975_deuterium_dv1104.par\",kFALSE,kFALSE,kFALSE,\"proton_D2_300torr.txt\",kTRUE,kTRUE,kTRUE,kFALSE,0,kTRUE,0.1,kFALSE,5.0,5.0,1,80.0,8)" > $O/logs/prod/$r.log 2>&1 < /dev/null
  vr=$(root -b -l -q "pid/root_ok.C(\"$f\")" 2>/dev/null | grep -E '^(VALID|INVALID)' | head -1); if [[ "$vr" == VALID* ]] && grep -q "dE/dx from CATIMA" $O/logs/prod/$r.log; then touch $f.done; echo "okrun $r $vr $(date +%T)"; else echo "FAIL $r ${vr:-unreadable}"; rm -f $f; fi; }
export -f one; export G
ls /mnt/f/a1975/gf_dp_find805/ | sed -nE 's/^(run_[0-9]+)_multifit_genfitter_p\.root$/\1/p' | xargs -P 6 -I{} bash -c 'one {}'
echo PROD_ALLDONE
