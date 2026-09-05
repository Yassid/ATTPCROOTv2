#!/usr/bin/env bash
# 46Ar(3He,d)47K: FIELD x PAD PLANE x STATE, in 3He + 5% CO2 at 300 torr.
#
#   ./run_3Hed_matrix.sh [-j N] [nEvents]
#
# Twelve samples: 2 fields (2.85, 3.9 T) x 2 pad planes (real AT-TPC, uniform 2 mm) x 3 states of
# 47K (1/2+ g.s., 0.36 MeV 3/2+, 2.02 MeV 7/2-). Pure 3He was on this matrix and was DROPPED on
# 2026-09-04 -- see the note at the end of this header.
#
# GENERATION IS SHARED, RECO IS NOT. Ar46_3Hed_sim.C takes the field as an argument and never
# reads the .par, so Geant4 transport depends on the FIELD and the GAS, not on the pad plane or on
# the transport numbers. Consequences, both exploited here:
#   * the two pad arms of one field re-reconstruct the SAME _sim.root, so a difference between
#     them is the pad plane and nothing else;
#   * the six 2.85 T generations already exist in /mnt/f/ar46_3hed_OLD_2.85T_placeholder (that
#     directory is "OLD" only in its placeholder-transport reco -- its sims are live), so only
#     3.9 T is generated from scratch. Do not move or rename it while this runs.
# Three generations instead of twelve is most of the cost of this campaign.
#
# TWO WAVES. Wave 1 is the AT-TPC pad plane at both fields, which is the adopted detector and
# therefore the headline: it answers "does 3.9 T beat 2.85 T for all three states" on its own.
# Wave 2 adds the 2 mm plane, which is the design question, and by then every sim it needs exists.
#
# THE 2 mm ARM CARRIES TWO KNOWN CONFOUNDS. Neither is fixed here; both are MEASURED, because a
# pitch comparison that silently measures something else is worse than none:
#   1. GAIN. A 2x2 mm pad collects 1/24 of the charge of an 8x12 mm one, so at a fixed PSA
#      threshold a fine plane simply loses hits. Measured on 14C(d,p): at gain 35000 the AT-TPC
#      plane gave 1870 usable protons against 848 for 2 mm; at 150000 they were 898 and 897. The
#      46Ar par is already at 150000, i.e. in the regime where the two planes are NOT separated by
#      the threshold -- but the per-arm hit and pattern-track counts are printed at the end so
#      that stays checkable rather than assumed.
#   2. PATTERN RECOGNITION. clusterRadius 15 mm and clusterDistance 7.5 mm were set for 8x12 mm
#      pads, where they span one to two pads; at 2 mm pitch they span seven. They are NOT re-tuned
#      (run_reco_Ar46_TC.C says so in its own header). So compare resolutions only after checking
#      that pattern efficiency did not move -- summarise_matrix prints it per arm.
#
# measSigma STAYS AT 0.6 mm IN BOTH ARMS even though the 2 mm plane would nominally want ~0.35.
# With no chi2 cut applied, measSigma is a pure renormalisation of chi2 and does not change the
# fitted trajectory (measured on 14C(d,p) at 7 T: chi2 scaled exactly as 1/measSigma^2 across
# 0.6/1.0/1.4 and the trajectory never moved). Holding it fixed keeps the pad plane the only
# difference between the arms; a per-arm value would confound the two.
#
# PURE 3He IS NOT HERE, and the reason is worth keeping. The readout window is hardware: AtPad
# fixes 512 time buckets (the AGET chip -- a par asking for 2048 is silently ignored) and
# AtDigiPar::GetTBTime accepts only 3/6/12/25/50/100 MHz, so the longest recordable drift is
# 512 tb @ 3 MHz = 163.84 us. Pure 3He ran v_d = 0.46-0.54 cm/us at 200 torr and 50-70 V/cm,
# i.e. 185-219 us for the metre -- it does not fit, and since pure 3He never saturates the fix is
# a much higher field cage, which is a proposal-level decision rather than a simulation setting.
set -uo pipefail
JOBS=4
if [ "${1:-}" = "-j" ]; then JOBS=${2:?-j needs a number}; shift 2; fi
NEV=${1:-12000}
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
ROOT_OUT=${MX_ROOT:-/mnt/f}
MASTER=$ROOT_OUT/ar46_3hed_matrix.log
say() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$MASTER"; }

# ONE DRIVER PER OUTPUT TREE. Two drivers sharing these directories do not merely duplicate work:
# each one's resume check reads the other's half-written output and takes it for a finished stage.
# The stale-lock branch tests the recorded PID, so a driver killed mid-run does not block forever.
DRVLOCK=$ROOT_OUT/ar46_3hed_matrix.lock
if ! mkdir "$DRVLOCK" 2>/dev/null; then
   other=$(cat "$DRVLOCK/pid" 2>/dev/null || echo "?")
   if [ "$other" != "?" ] && kill -0 "$other" 2>/dev/null; then
      echo "REFUSING TO START: matrix driver PID $other is already running"; exit 2
   fi
   echo "removing stale lock from PID ${other:-?}"; rm -rf "$DRVLOCK"; mkdir "$DRVLOCK" || exit 2
fi
echo $$ > "$DRVLOCK/pid"
trap 'rm -rf "$DRVLOCK"' EXIT INT TERM

# field : par suffix : sim directory (generation lives here, shared by both pad arms)
FIELDS=(
  "2.85:B285:/mnt/f/ar46_3hed_OLD_2.85T_placeholder"
  "3.9:B39:/mnt/f/ar46_3hed_gen_B39"
)
# pad tag : padSize_mm argument (<=0 keeps the real AT-TPC plane)
PADS=( "attpc:-1" "2mm:2.0" )
# state : seed -- distinct seeds per state, or two states share a random sequence and the
# "independent" samples are the same events with a different Q value.
STATES=( "gs:3001" "360:3011" "2020:3021" )

for f in "${FIELDS[@]}"; do
  IFS=: read -r BT SFX SIMDIR <<< "$f"
  P="$REPO/parameters/ATTPC.46Ar_3Hed_sim_${SFX}.par"
  [ -f "$P" ] || { echo "MISSING $P -- run make_3Hed_pars.sh (needs the Magboltz point for $BT T)"; exit 2; }
  grep -q "MAGBOLTZ" "$P" || { echo "$P carries no Magboltz transport"; exit 2; }
  pb=$(awk '/^BField:Double_t/{print $2}' "$P")
  awk -v a="$pb" -v b="$BT" 'BEGIN{exit !(a-b<1e-6 && b-a<1e-6)}' \
     || { echo "$P says BField $pb but this arm is $BT T"; exit 2; }
done

# --- one sample -------------------------------------------------------------------------------
# accumulate_3Hed.sh already does generate -> reco -> pid, is resumable per stage on evidence, and
# takes the pad size and a separate sim directory. The fit is here because accumulate deliberately
# stops before it (the PID gate is drawn by hand on the plane it produces).
one() { # field parSuffix simDir padTag padArg state seed
  local BT=$1 SFX=$2 SIMDIR=$3 PTAG=$4 PARG=$5 ST=$6 SEED=$7
  local OUT=$ROOT_OUT/ar46_3hed_mx_${SFX}_${PTAG}
  local J=${ST}_s${SEED}
  mkdir -p "$OUT"
  cat > "$OUT/PROVENANCE.txt" <<EOF
FIELD x PAD MATRIX arm ${SFX}_${PTAG}. 46Ar(3He,d)47K in 3He+5% CO2, 300 torr, 293.15 K.
B = $BT T, pad plane = $PTAG ($PARG mm; <=0 means the real AT-TPC plane with the beam hole).
Operating point E = 140 V/cm, 512 tb @ 6 MHz, per-field Magboltz transport from
parameters/ATTPC.46Ar_3Hed_sim_${SFX}.par. genfit at NEGATIVE bField, measSigma 0.6 (the SAME
value in both pad arms, on purpose), matEffects + CATIMA, backward seed fix, back-extrapolation.
Generation shared with the other pad arm of this field, in $SIMDIR.
Produced by run_3Hed_matrix.sh. Map: /mnt/f/AR46_3HED_INDEX.md
EOF

  if [ ! -f "$OUT/$J.marker" ]; then
    say "$SFX/$PTAG $J accumulate (B = $BT T, pads $PARG)"
    "$HERE/accumulate_3Hed.sh" "$ST" "$SEED" "$NEV" "$OUT" "$BT" "$PARG" "$SIMDIR" \
       >> "$MASTER" 2>&1 || { say "$SFX/$PTAG $J ACCUMULATE_FAILED"; return 1; }
  fi

  if grep -q "Done" "$OUT/${J}_fit.log" 2>/dev/null && [ -s "$OUT/${J}_genfitter_d.root" ]; then
    say "$SFX/$PTAG $J already fitted"
  else
    say "$SFX/$PTAG $J genfit"
    ( set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
      export VMCWORKDIR="$REPO"
      export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
      cd "$HERE"
      root -b -q -l "fitGenfitter_Ar46.C(\"$J\",-1,\"$OUT/\",\"\",\"\",-${BT},2,5,\"\",0.6,10.0,170.0,kTRUE,kTRUE,\"ATTPC_He3CO2_300torr\",\"ATTPC.46Ar_3Hed_sim_${SFX}.par\")" \
        > "$OUT/${J}_fit.log" 2>&1 ) || true
    grep -q "Done" "$OUT/${J}_fit.log"           || { say "$SFX/$PTAG $J FIT_FAILED"; return 1; }
    grep -q "dE/dx from CATIMA" "$OUT/${J}_fit.log" || { say "$SFX/$PTAG $J CATIMA_NOT_ENABLED"; return 1; }
  fi

  # A fit that produced almost nothing is a configuration failure, not a physics result -- say so
  # here rather than letting it reach the summary as a suspiciously good resolution on 12 tracks.
  local nfit
  nfit=$(grep -c "tracks fitted" "$OUT/${J}_fit.log" 2>/dev/null || echo 0)
  say "$SFX/$PTAG $J DONE ($nfit events with a fitted track)"
  [ "$nfit" -ge 500 ] || say "$SFX/$PTAG $J *** ONLY $nfit FITTED EVENTS -- check this arm before using it ***"
}

# --- waves ------------------------------------------------------------------------------------
wave() { # padTag padArg
  local PTAG=$1 PARG=$2
  local pids=()
  for f in "${FIELDS[@]}"; do
    IFS=: read -r BT SFX SIMDIR <<< "$f"
    mkdir -p "$SIMDIR"
    for s in "${STATES[@]}"; do
      IFS=: read -r ST SEED <<< "$s"
      while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n; done
      one "$BT" "$SFX" "$SIMDIR" "$PTAG" "$PARG" "$ST" "$SEED" &
      pids+=($!)
      sleep 15   # stagger the geometry open, the one contended step
    done
  done
  local bad=0
  for p in "${pids[@]}"; do wait "$p" || bad=1; done
  return $bad
}

say "########## matrix start: ${#FIELDS[@]} fields x ${#PADS[@]} planes x ${#STATES[@]} states, $NEV entries, -j$JOBS ##########"
fail=0
for p in "${PADS[@]}"; do
  IFS=: read -r PTAG PARG <<< "$p"
  say "===== WAVE: $PTAG pad plane ====="
  wave "$PTAG" "$PARG" || { fail=1; say "wave $PTAG had failures"; }
done
say "########## matrix complete (fail=$fail) ##########"

cat <<'EOF'

Summarise the field axis at each pad plane with:

  root -b -q 'ex_genfit_3Hed.C("gs_s3001,360_s3011,2020_s3021",
     "/mnt/f/ar46_3hed_mx_B285_attpc,/mnt/f/ar46_3hed_mx_B39_attpc",
     "2.85 T AT-TPC,3.9 T AT-TPC",10.0,100.0,kTRUE,-1.0,
     "/mnt/f/ar46_3hed_OLD_2.85T_placeholder,/mnt/f/ar46_3hed_gen_B39")'

and the pad axis at each field by swapping the _attpc arms for _2mm. NOTE that ex_genfit_3Hed.C
pools every tag it is given into one distribution, so run it ONE STATE AT A TIME when you want a
per-state resolution -- three states pooled is a three-peak spectrum, and its IQR is a statement
about the level spacing, not about the detector.
EOF
exit $fail
