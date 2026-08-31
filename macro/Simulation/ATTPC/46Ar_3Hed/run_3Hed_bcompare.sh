#!/usr/bin/env bash
# 46Ar(3He,d)47K GROUND STATE at 2.0 / 2.85 / 3.8 T, on the real AT-TPC pad plane WITH the beam
# hole (a telescope goes behind it for the heavy fragment, so the hole is not negotiable).
#
#   ./run_3Hed_bcompare.sh
#
# SINGLE CORE THROUGHOUT -- the 14C(d,p) campaign owns the rest of the machine.
#
# WHAT IS NEW HERE, and why the earlier field comparison does not answer this question:
#   * Every field gets its OWN transport. DriftVelocity and CoefL/CoefT now come from a Magboltz
#     scan of 3He+5%CO2 at 300 torr, 50 V/cm, AT THAT FIELD (make_3Hed_pars.sh). The August
#     comparison held them at a1975 H2 placeholders identical across fields, so it varied only
#     the magnetic curvature and understated the higher field, whose transverse diffusion really
#     is suppressed.
#   * The genfit fit is the post-54a3726c one, at NEGATIVE bField -- see fitGenfitter_Ar46.C.
#
# GENERATION IS REUSED WHERE IT EXISTS AND IS STILL VALID. Ar46_3Hed_sim.C takes the field as an
# argument and never reads the .par, so the Geant4 transport at 2.85 and 3.8 T is unaffected by
# the new transport numbers -- only digitisation is. Those two fields therefore re-reconstruct
# from their existing _sim.root files and only 2.0 T is generated from scratch.
#
# OUTPUT GOES TO FRESH _mb_ DIRECTORIES so Magboltz-transport products can never be confused with
# the placeholder-transport ones still on disk. Do not point these at the old directories.
set -eo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
NEV=${NEV:-12000}
SEEDS=${SEEDS:-"3001 3002"}
MASTER=/mnt/f/ar46_3hed_bcompare.log

say() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$MASTER"; }

# field : output dir : sim dir (where _sim.root lives / goes) : par suffix
CFG=(
  "2.0:/mnt/f/ar46_3hed_mb_B20:/mnt/f/ar46_3hed_mb_B20:B20"
  "2.85:/mnt/f/ar46_3hed_mb_B285:/mnt/f/ar46_3hed:B285"
  "3.8:/mnt/f/ar46_3hed_mb_B38:/mnt/f/ar46_3hed_B38:B38"
)

# Refuse to start unless every par exists AND actually carries Magboltz numbers. A par still
# holding the 1.30 placeholder would run silently and produce a comparison that looks fine.
for c in "${CFG[@]}"; do
  IFS=: read -r BT OUT SIMDIR SFX <<< "$c"
  P="$REPO/parameters/ATTPC.46Ar_3Hed_sim_${SFX}.par"
  [ -f "$P" ] || { echo "MISSING $P -- run make_3Hed_pars.sh"; exit 2; }
  grep -q "MAGBOLTZ" "$P" || { echo "$P carries no Magboltz DriftVelocity -- re-run make_3Hed_pars.sh"; exit 2; }
  # and the field in the par must be the field we are about to simulate
  pb=$(awk '/^BField:Double_t/{print $2}' "$P")
  awk -v a="$pb" -v b="$BT" 'BEGIN{exit !(a-b<1e-6 && b-a<1e-6)}' \
     || { echo "$P says BField $pb but this arm is $BT T"; exit 2; }
done
say "########## g.s. field comparison start: 2.0 / 2.85 / 3.8 T, $NEV entries, seeds $SEEDS ##########"

for c in "${CFG[@]}"; do
  IFS=: read -r BT OUT SIMDIR SFX <<< "$c"
  mkdir -p "$OUT" "$SIMDIR"
  for s in $SEEDS; do
    t0=$SECONDS
    say "--- B = $BT T, seed $s (sims $SIMDIR -> $OUT) ---"
    "$DIR/accumulate_3Hed.sh" gs "$s" "$NEV" "$OUT" "$BT" -1 "$SIMDIR" 2>&1 | tee -a "$MASTER"
    say "--- B = $BT T seed $s accumulated in $(( (SECONDS-t0)/60 )) min ---"
    # fit immediately, at NEGATIVE field and this arm's own par
    J=gs_s${s}
    if grep -q "Done" "$OUT/${J}_fit.log" 2>/dev/null && [ -s "$OUT/${J}_genfitter_d.root" ]; then
      say "$J already fitted"
    else
      say "$J genfit"
      ( cd "$REPO" && export VMCWORKDIR=$PWD && set +u && source build/config.sh >/dev/null 2>&1 && set -u
        root -b -q -l "macro/Simulation/ATTPC/46Ar_3Hed/fitGenfitter_Ar46.C(\"$J\",-1,\"$OUT/\",\"\",\"\",-${BT},2,5,\"\",0.6,10.0,170.0,kTRUE,kTRUE,\"ATTPC_He3CO2_300torr\",\"ATTPC.46Ar_3Hed_sim_${SFX}.par\")" \
          > "$OUT/${J}_fit.log" 2>&1 )
      grep -q "Done" "$OUT/${J}_fit.log" || { say "$J FIT_FAILED"; exit 1; }
      grep -q "dE/dx from CATIMA" "$OUT/${J}_fit.log" || { say "$J CATIMA_NOT_ENABLED"; exit 1; }
    fi
  done
done
say "########## field comparison complete ##########"
