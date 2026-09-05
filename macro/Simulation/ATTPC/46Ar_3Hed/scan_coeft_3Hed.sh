#!/usr/bin/env bash
# TRANSVERSE-DIFFUSION SCAN at the adopted 2.85 T operating point.
#
#   ./scan_coeft_3Hed.sh
#
# WHY THIS EXISTS. Going from the placeholder transport to the real Magboltz one made sigma(Ex)
# WORSE, not better: 0.258 -> 0.415 MeV at theta_lab 125-145 and 0.465 -> 0.594 at 110-125. Less
# diffusion should not cost resolution, so either the mechanism is pad under-sampling -- real D_T
# is 1.25 mm rms at full drift against 8 x 12 mm pads, so the charge lands on too few pads for
# charge-weighted centroiding to beat the pitch, and the placeholder's excess diffusion was
# accidentally helping -- or the effect is not diffusion at all.
#
# THE PLACEHOLDER COMPARISON CANNOT SETTLE IT, because those two pars differ in FOUR lines:
# EField 5000 -> 14000 V/m, DriftVelocity 1.30 -> 1.3196, CoefL 0.0009 -> 0.00060951 and
# CoefT 0.0009 -> 0.000103307. This scan changes ONE: CoefT, from the adopted B285 par.
# CoefL is deliberately held fixed -- the hypothesis is transverse under-sampling of the pad
# plane, and moving both would not say which one acted.
#
# NO NEW GENERATION. Geant4 transport depends on the B field, never on the .par, so every point
# re-digitises and re-reconstructs the SAME _sim.root that the adopted 2.85 T arm used. The
# ct10 point is not run here: it IS /mnt/f/ar46_3hed_mb_B285, bit-for-bit the same configuration.
#
# THE SCAN VARIABLE IS D_T, NOT CoefT. AtClusterize uses sigma_T[mm] = 10*sqrt(2*CoefT*t_drift),
# so CoefT goes as D_T^2 and a list of round CoefT values would be a lopsided list of widths.
# Multipliers below are on D_T; sigma_T is quoted at the full 100 cm drift (75.8 us at v_d 1.3196).
#
#   tag    xD_T   CoefT [cm^2/us]   D_T [um/sqrt(cm)]   sigma_T at 100 cm [mm]
#   ct05   0.5    2.582675e-05       62.6                0.63
#   ct10   1.0    1.033070e-04      125.1                1.25   <- ADOPTED, already on disk
#   ct20   2.0    4.132280e-04      250.2                2.50
#   ct30   3.0    9.297630e-04      375.3                3.75   (~ the old placeholder, 0.0009)
#   ct40   4.0    1.652912e-03      500.4                5.00
#
# Resumable per stage on EVIDENCE, not on file existence: a written file is not a finished job.
set -eo pipefail
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO=/home/yassid/fair_install/ATTPCROOTv2-OpenKF
SIMDIR=/mnt/f/ar46_3hed_OLD_2.85T_placeholder     # holds the live 2.85 T _sim.root
BASEPAR="$REPO/parameters/ATTPC.46Ar_3Hed_sim_B285.par"
TAG=${TAG:-gs_s3001}
BT=2.85
JOBS=${JOBS:-4}
MASTER=/mnt/f/ar46_3hed_coeft.log

say() { echo "[$(date +%F' '%H:%M:%S)] $*" | tee -a "$MASTER"; }

# tag : CoefT
POINTS=(
  "ct05:2.582675e-05"
  "ct20:4.132280e-04"
  "ct30:9.297630e-04"
  "ct40:1.652912e-03"
)

[ -f "$BASEPAR" ] || { echo "MISSING $BASEPAR -- run make_3Hed_pars.sh"; exit 2; }
grep -q "MAGBOLTZ" "$BASEPAR" || { echo "$BASEPAR is not the Magboltz par"; exit 2; }
[ -s "$SIMDIR/${TAG}_sim.root" ] || { echo "MISSING $SIMDIR/${TAG}_sim.root"; exit 2; }
# the base par must be the one the reference point was built with, or the scan has no anchor
base_ct=$(awk '/^CoefT:Double_t/{print $2}' "$BASEPAR")
awk -v a="$base_ct" 'BEGIN{exit !(a-1.03307e-4<1e-9 && 1.03307e-4-a<1e-9)}' \
   || { echo "base par CoefT is $base_ct, expected 1.03307e-04 -- the scan anchor moved"; exit 2; }

# --- pars: the B285 par with CoefT replaced, and nothing else touched -------------------------
for p in "${POINTS[@]}"; do
  IFS=: read -r t ct <<< "$p"
  P="$REPO/parameters/ATTPC.46Ar_3Hed_sim_B285_${t}.par"
  awk -v ct="$ct" -v tag="$t" '
    /^CoefT:Double_t/ { printf "CoefT:Double_t %20s   # cm^2/us  DIFFUSION SCAN point %s -- CoefT is the ONLY line changed from ATTPC.46Ar_3Hed_sim_B285.par\n", ct, tag; next }
    { print }
  ' "$BASEPAR" > "$P"
  got=$(awk '/^CoefT:Double_t/{print $2}' "$P")
  [ "$got" = "$ct" ] || { echo "$P did not take CoefT $ct (got $got)"; exit 2; }
  # nothing but that one line may differ
  nd=$(diff <(grep -v '^CoefT:Double_t' "$BASEPAR") <(grep -v '^CoefT:Double_t' "$P") | wc -l)
  [ "$nd" -eq 0 ] || { echo "$P differs from the base par in more than CoefT"; exit 2; }
done
say "########## CoefT scan start: $TAG, ${#POINTS[@]} points, $JOBS at a time ##########"

run_point() {
  local t=$1 ct=$2
  local OUT=/mnt/f/ar46_3hed_${t}
  local PAR=ATTPC.46Ar_3Hed_sim_B285_${t}.par
  mkdir -p "$OUT"
  set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
  export VMCWORKDIR="$REPO"
  export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
  cd "$DIR"

  cat > "$OUT/PROVENANCE.txt" <<EOF
DIFFUSION SCAN point $t. 46Ar(3He,d)47K ground state, 2.85 T, adopted operating point
(E = 140 V/cm, v_d 1.3196 cm/us, 512 tb @ 6 MHz, AT-TPC pad plane with the beam hole).
CoefT = $ct cm^2/us -- THE ONLY LINE that differs from ATTPC.46Ar_3Hed_sim_B285.par.
Reconstructed from $SIMDIR/${TAG}_sim.root, the same sim the adopted 2.85 T arm used.
The x1 point of this scan is /mnt/f/ar46_3hed_mb_B285 itself; it is not duplicated here.
Produced by scan_coeft_3Hed.sh. Map: /mnt/f/AR46_3HED_INDEX.md
EOF

  # --- reco (this is where CoefT acts: AtClusterize sigma_T = 10*sqrt(2*CoefT*t)) -------------
  if [ -s "$OUT/${TAG}_reco.root" ] && grep -q "sim reco done" "$OUT/${TAG}_reco.log" 2>/dev/null; then
    say "$t reco already complete"
  else
    say "$t reco (CoefT $ct)"
    root -b -q -l "run_reco_Ar46_TC.C(\"$SIMDIR/${TAG}_sim.root\",\"$OUT/${TAG}_reco.root\",\"$PAR\",20,0,20.0,15.0,7.5,\"tc\",20,-1)" \
         > "$OUT/${TAG}_reco.log" 2>&1 || true
    grep -q "sim reco done" "$OUT/${TAG}_reco.log" || { say "$t RECO_FAILED"; return 1; }
    # the par really reached the clusterizer -- AtClusterize prints what it read
    grep -q "Transverse coefficient of diffusion: $(awk -v c="$ct" 'BEGIN{printf "%g", c}')" "$OUT/${TAG}_reco.log" \
      || say "$t WARNING: could not confirm CoefT in the reco log -- check it by hand"
  fi

  # --- fit, identical in every respect to the adopted arm except the par ----------------------
  if grep -q "Done" "$OUT/${TAG}_fit.log" 2>/dev/null && [ -s "$OUT/${TAG}_genfitter_d.root" ]; then
    say "$t fit already complete"
  else
    say "$t genfit"
    root -b -q -l "fitGenfitter_Ar46.C(\"$TAG\",-1,\"$OUT/\",\"\",\"\",-${BT},2,5,\"\",0.6,10.0,170.0,kTRUE,kTRUE,\"ATTPC_He3CO2_300torr\",\"$PAR\")" \
         > "$OUT/${TAG}_fit.log" 2>&1 || true
    grep -q "Done" "$OUT/${TAG}_fit.log" || { say "$t FIT_FAILED"; return 1; }
    grep -q "dE/dx from CATIMA" "$OUT/${TAG}_fit.log" || { say "$t CATIMA_NOT_ENABLED"; return 1; }
  fi
  echo COMPLETED > "$OUT/${TAG}.marker"
  say "$t DONE"
}

pids=()
for p in "${POINTS[@]}"; do
  IFS=: read -r t ct <<< "$p"
  while [ "$(jobs -rp | wc -l)" -ge "$JOBS" ]; do wait -n; done
  run_point "$t" "$ct" &
  pids+=($!)
  sleep 20   # stagger: four ROOT jobs opening the same geometry at once is the one contended step
done
fail=0
for pid in "${pids[@]}"; do wait "$pid" || fail=1; done
say "########## CoefT scan complete (fail=$fail) ##########"

O=/mnt/f/ar46_3hed_OLD_2.85T_placeholder
cat <<EOF

Analyse with (the x1 point is the adopted arm; every point re-reconstructs the SAME sim, so the
sim directory is repeated once per arm -- ex_genfit_3Hed.C falls back to the fit's own directory
for any arm the list does not reach, and there is no _sim.root there):

  O=/mnt/f/ar46_3hed_OLD_2.85T_placeholder
  root -b -q "ex_genfit_3Hed.C(\"gs_s3001\",
     \"/mnt/f/ar46_3hed_ct05,/mnt/f/ar46_3hed_mb_B285,/mnt/f/ar46_3hed_ct20,/mnt/f/ar46_3hed_ct30,/mnt/f/ar46_3hed_ct40\",
     \"0.5x D_T,1x (adopted),2x,3x,4x\",10.0,100.0,kTRUE,-1.0,
     \"$O,$O,$O,$O,$O\")"
EOF
exit $fail
