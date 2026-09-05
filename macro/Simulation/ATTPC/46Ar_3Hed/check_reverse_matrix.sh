#!/usr/bin/env bash
# Two checks on the ReverseDrift feature that the first validation did not cover.
#
#   ./check_reverse_matrix.sh          # both
#   PART=1 ./check_reverse_matrix.sh   # regression only
#   PART=2 ./check_reverse_matrix.sh   # matrix only
#
# PART 1 -- NORMAL MODE IS UNCHANGED, MEASURED RATHER THAN ARGUED.
# The claim so far was structural: CalculateZGeo early-returns before the reversal and
# AtClusterize picks the original formula with a ternary, so nothing in normal mode takes a new
# code path. That is a code-reading argument, and code-reading arguments are how silent breakage
# gets in. This re-runs a full production reconstruction with the CURRENT build and diffs it
# against the product already on disk, hit by hit.
#
# This is only a meaningful diff because THE DIGITISATION IS REPRODUCIBLE: nothing in the chain
# calls SetSeed, so gRandom keeps ROOT's default TRandom3 sequence and two runs of the same input
# draw the same numbers. If a seed is ever introduced, this check has to become statistical and
# this comment is the warning that it silently would not be a regression test any more.
#
# PART 2 -- THE MAPPING HOLDS ACROSS THE MATRIX, not just where it was first tested.
# The truth-anchored validation was run on ONE cell: 2.85 T, AT-TPC pads, 400 events. The pad
# pitch changes the pad response and the field changes the transverse diffusion, and both feed the
# digitisation, so the invariant is re-checked in every cell rather than assumed to transfer.
#
# Note on PSA coverage: AtPSAMax already reaches AtPSA::CalibrateZ (AtPSAMax.cxx:50), so the
# calibration entry point IS exercised here. The paths still untested are the seven PSAs that call
# CalculateZGeo directly -- covered by construction, since that is the function the flip lives in --
# and CalibrateZ's Spyral branch, which no production macro enables and which is deliberately left
# on its own (pre-existing, opposite) z convention. See docs/REVERSED_DETECTOR.md.
set -uo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
OUT=${OUT:-/mnt/f/ar46_3hed_revcheck}
NEV2=${NEV2:-400}
PART=${PART:-0}
mkdir -p "$OUT"
LOG="$OUT/check.log"
say() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$LOG"; }

set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export VMCWORKDIR="$REPO"
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
cd "$HERE"

A=/mnt/f/ar46_3hed_OLD_2.85T_placeholder
B=/mnt/f/ar46_3hed_gen_B39

# ---------------------------------------------------------------- PART 1
if [ "$PART" = "0" ] || [ "$PART" = "1" ]; then
  say "########## PART 1: normal-mode regression against the on-disk production ##########"
  REF=/mnt/f/ar46_3hed_mx_B285_attpc/gs_s3001_reco.root
  NEW=$OUT/regress_normal_reco.root
  if [ ! -s "$REF" ]; then
    say "PART 1 SKIPPED: reference $REF is missing"
  else
    if [ -s "$NEW" ]; then
      say "regression reco already built"
    else
      say "re-running the 2.85 T AT-TPC production reco with the current build (12000 entries, ~80 min)"
      root -b -q -l "run_reco_Ar46_TC.C(\"$A/gs_s3001_sim.root\",\"$NEW\",\"ATTPC.46Ar_3Hed_sim_B285.par\",20,0,20.0,15.0,7.5,\"tc\",20,-1)" \
           > "$OUT/regress_normal.log" 2>&1 || true
      grep -q "sim reco done" "$OUT/regress_normal.log" || { say "PART 1 RECO FAILED -- not comparing"; rm -f "$NEW"; }
    fi
    # Test the PRODUCT, not the filename: a killed reco leaves a large, non-empty, UNCLOSED file,
    # and handing that to the comparison crashed it rather than reporting anything useful.
    if [ -s "$NEW" ] && grep -q "sim reco done" "$OUT/regress_normal.log"; then
      say "diffing hit-by-hit against $REF"
      root -b -q -l "compare_reco_3Hed.C(\"$REF\",\"$NEW\")" 2>&1 | tee -a "$LOG" | grep -E "IDENTICAL|DIFFER|events|hits"
    fi
  fi
fi

# ---------------------------------------------------------------- PART 2
if [ "$PART" = "0" ] || [ "$PART" = "2" ]; then
  say "########## PART 2: truth-anchored validation across field x pad plane ##########"
  # tag : sim dir : par stem : pad size (mm, <=0 = real AT-TPC plane)
  CFG=(
    "b285_attpc:$A:ATTPC.46Ar_3Hed_sim_B285:-1"
    "b285_2mm:$A:ATTPC.46Ar_3Hed_sim_B285:2.0"
    "b39_attpc:$B:ATTPC.46Ar_3Hed_sim_B39:-1"
    "b39_2mm:$B:ATTPC.46Ar_3Hed_sim_B39:2.0"
  )
  for c in "${CFG[@]}"; do
    IFS=: read -r tag sim stem pad <<< "$c"
    [ -s "$sim/gs_s3001_sim.root" ] || { say "$tag SKIPPED (no sim in $sim)"; continue; }
    for m in nrm rev; do
      P="${stem}.par"; [ "$m" = rev ] && P="${stem}_rev.par"
      [ -f "$REPO/parameters/$P" ] || { say "$tag/$m SKIPPED (no par $P)"; continue; }
      f="$OUT/${tag}_${m}_reco.root"
      # Resume on EVIDENCE, not on the filename: a killed job leaves a large unclosed file that
      # every size test passes. Same trap as PART 1.
      if [ -s "$f" ] && grep -q "sim reco done" "$OUT/${tag}_${m}.log" 2>/dev/null; then
        say "$tag/$m reco already built"; continue
      fi
      rm -f "$f"
      say "$tag/$m reco ($NEV2 entries, pads $pad, par $P)"
      root -b -q -l "run_reco_Ar46_TC.C(\"$sim/gs_s3001_sim.root\",\"$f\",\"$P\",20,$NEV2,20.0,15.0,7.5,\"tc\",20,$pad)" \
           > "$OUT/${tag}_${m}.log" 2>&1 || true
      grep -q "sim reco done" "$OUT/${tag}_${m}.log" || say "$tag/$m RECO FAILED"
    done
    if grep -q "sim reco done" "$OUT/${tag}_nrm.log" 2>/dev/null &&
       grep -q "sim reco done" "$OUT/${tag}_rev.log" 2>/dev/null; then
      say "--- validating $tag ---"
      root -b -q -l "validate_reverse_3Hed.C(\"$sim/gs_s3001_sim.root\",\"$OUT/${tag}_nrm_reco.root\",\"$OUT/${tag}_rev_reco.root\")" \
        2>&1 | tee -a "$LOG" | grep -E "PASS|FAIL|median|normal|reversed|control"
    fi
  done
fi
say "########## checks complete ##########"
