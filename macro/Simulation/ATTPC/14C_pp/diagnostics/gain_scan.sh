#!/usr/bin/env bash
# Gain scan on CLEAN protons. Targets, measured on a1954 data (gated protons, production reco):
#     charge/mm ~1113      pads/mm ~1.61      hits-per-pad ~1.02
# Current sim at Gain 150000 gives charge/mm 65.6 and pads/mm 0.40.
#
# Measured on truth-gated protons ONLY. The earlier "gain is right" conclusion came from a sample
# where 67 % of proton clusters were merged with the beam, and 14C ionises ~36x more per mm than
# a proton, so the beam supplied nearly all the charge. With the 3 cm hole MERGED is now 0.1 %.
#
# Same fixed MC truth for every point (diagnostics/negB/attpcsim_negB.root), so only digitization
# changes. 1500 events is plenty for medians (~700 protons).
set -eo pipefail
G=$1
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd); SIM=$(cd "$HERE/.." && pwd)
REPO=$(cd "$SIM/../../../.." && pwd); OUT=/mnt/f/a1954_C14_gain; mkdir -p "$OUT"
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
PAR="ATTPC.a1954_C14_sim_g${G}.par"
# only ever REPLACE a value on an existing line: a standalone comment makes FairParAsciiFileIo
# silently drop every parameter below it
sed -e "s/^Gain:Double_t.*/Gain:Double_t               $G   # gain scan/" \
    "$REPO/parameters/ATTPC.a1954_C14_sim.par" > "$REPO/parameters/$PAR"
grep -qE "^Gain:Double_t +$G" "$REPO/parameters/$PAR" || { echo "g$G PAR_FAILED"; exit 1; }
cd "$SIM"
root -b -q -l "run_reco_C14.C(\"$HERE/negB/attpcsim_negB.root\",\"$OUT/g${G}_reco.root\",\"$PAR\",20,20,8,1500,30.0,\"mover\")" > "$OUT/g${G}_reco.log" 2>&1
[ -s "$OUT/g${G}_reco.root" ] || { echo "g$G RECO_FAILED"; exit 1; }
root -b -q -l "gate_truth_C14.C(\"$OUT/g${G}_reco.root\",\"$OUT/\",\"g${G}g\",0.9,10)" > "$OUT/g${G}_gate.log" 2>&1
[ -s "$OUT/g${G}g_reco.root" ] || { echo "g$G GATE_FAILED"; exit 1; }
"$HERE/match_diag" "$OUT/g${G}g_reco.root" "g$G" 1.30 -1 > "$OUT/g${G}_diag.txt" 2>&1
grep -q DECOMPOSE "$OUT/g${G}_diag.txt" || { echo "g$G DIAG_FAILED"; exit 1; }
rm -f "$OUT/g${G}_reco.root" "$OUT/g${G}g_reco.root"
echo "GAIN $G : $(grep LENGTH "$OUT/g${G}_diag.txt" | sed 's/.*hits\/mm/hits\/mm/') | $(grep DECOMPOSE "$OUT/g${G}_diag.txt" | cut -d: -f2)"
echo COMPLETED > "$OUT/g${G}.marker"
