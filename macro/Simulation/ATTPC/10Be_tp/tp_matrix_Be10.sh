#!/usr/bin/env bash
# The deliverable table: for every configuration and every theta_cm slice, the fitted core
# resolution, the significance of the 0+_2 against BOTH nulls, and whether its area ratio comes
# back at the value that was injected.
#
#   ./tp_matrix_Be10.sh [root_dir] [nTot]
#
# A-vs-B is the plain null (no 0+_2 at all). A-vs-C is the one that matters: no 0+_2, but the
# 2.109 free to slide and broaden -- the freedom an unresolved doublet would exploit. Quote C.
#
# The ratio column is the fitted 0+_2/2+ area against the ~0.19 actually injected (1/5 of the
# population, times the small difference in the two levels' acceptances). A significance without
# a recovered ratio is not a measurement of the state, it is a measurement of something.
set -uo pipefail
ROOTDIR=${1:-/mnt/f/Be10_tp}
NTOT=${2:-20000}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
cd "$HERE"; mkdir -p plots
SLICES="2:12 12:22 22:32 32:45 2:45"
CFGS="b285_attpc b285_2mm b400_attpc b400_2mm b700_attpc b700_2mm"

run() { # cfg lo hi vtx
   root -b -q -l "tp_spectrum_Be10.C(\"$ROOTDIR\",\"$1\",0.2,$2,$3,$NTOT,$4,112.20,\"plots\")" 2>&1 |
      sed 's/\x1b\[[0-9;]*m//g'
}

echo "############ 10Be(t,p)12Be : the 0+_2 at 2.251 MeV, 5x suppressed ############"
echo "root $ROOTDIR, $NTOT counts per slice, Asimov (expected) significances, sqrt(N) scaling"
for VTX in kTRUE kFALSE; do
   echo
   echo "===== beam energy: $([ $VTX = kTRUE ] && echo 'AT THE RECONSTRUCTED VERTEX' || echo 'CONSTANT') ====="
   printf "%-12s %-8s %9s %9s %9s %10s %8s\n" config slice sig_core "A-vs-B" "A-vs-C" ratio verdict
   for cfg in $CFGS; do
      ls "$ROOTDIR/$cfg"/exres_ex2251_s*_"$cfg".root >/dev/null 2>&1 || { printf "%-12s %s\n" "$cfg" "(not yet run)"; continue; }
      for sl in $SLICES; do
         lo=${sl%%:*} hi=${sl##*:}
         out=$(run "$cfg" "$lo" "$hi" "$VTX")
         sc=$(echo "$out" | grep -oP 'core sigma \K[0-9.]+')
         ab=$(echo "$out" | grep -oP 'expected significance \K[0-9.]+')
         ac=$(echo "$out" | grep -oP 'A beats it by [-0-9.]+  =  \K[0-9.]+')
         rt=$(echo "$out" | grep -oP 'area ratio \K[0-9.]+')
         # "recovered"/"NOT recovered" comes from the macro, which compares against the ratio it
         # actually injected in THAT slice rather than against a nominal 0.2
         vd=$(echo "$out" | grep -oP '\(\K(NOT recovered|recovered)')
         [ -z "${vd:-}" ] && vd=$(echo "$out" | grep -q "reference peak collapsed" && echo "COLLAPSED" || echo "-")
         printf "%-12s %-8s %9s %9s %9s %10s %8s\n" "$cfg" "$lo-$hi" "${sc:--}" "${ab:--}" "${ac:--}" "${rt:--}" "$vd"
      done
   done
done
echo
echo "matrix done"
