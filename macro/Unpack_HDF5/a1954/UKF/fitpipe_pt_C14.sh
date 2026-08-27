#!/usr/bin/env bash
# a1954 14C(p,t)12C: gate on the triton, fit with GENFIT+CATIMA, cache the kinematics.
#
# WHY THIS IS NOT THE (p,p') PIPELINE WITH A DIFFERENT GATE:
#   * the ejectile is a triton, so the fitter needs particle="triton" -- pdg 1000010030 and the
#     triton mass. Fitting a triton as a proton gets the rigidity right and the ENERGY wrong,
#     which is exactly the quantity the excitation energy is built from.
#   * the residual is 12C, not 14C, and Q = -4.641 MeV. At Ecm = 10.72 MeV that leaves only
#     6.08 MeV of excitation available, so 12C states above that are CLOSED -- a useful check:
#     the 7.654 and 9.641 loci should be empty, and on the gated PID plane they are.
#   * the triton emerges at theta_lab 8-24 deg with Brho ~ 1.5 Tm, a helix radius of ~52 cm in a
#     29 cm chamber, so these tracks LEAVE rather than curl and their rigidity comes from a partial
#     arc. Expect worse resolution than the proton channel, and check it rather than assuming it.
#
#   ./fitpipe_pt_C14.sh "run_0055 run_0058 run_0060" 3
set -eo pipefail
RUNS="${1:-run_0055}"; NPAR="${2:-3}"
HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
REPO=$(cd "$HERE/../../../.." && pwd)
IN="${IN:-/mnt/f/a1954_C14_reco_hdb_slim/}"
GATEDIR="${GATEDIR:-/mnt/f/a1954_C14_pt_in/}"
OUTDIR="${OUTDIR:-/mnt/f/a1954_C14_pt_fit/}"
# POLAR WINDOW GIVEN TO THE FITTER, and it is not the (p,p') one.
# fitGenfit_C14.C applies SetThetaWindow(THMIN, THMAX) on the RECONSTRUCTED polar, and the
# reconstruction reports 180 - theta_lab. The proton analysis uses [10,170], which for a recoil
# proton (reconstructed 90-180) only clips theta_cm > 160 and is harmless. The triton band is
# reconstructed 156-172, so an upper edge at 170 SLICES THROUGH IT: 998 of 8121 band tracks, 12 %,
# and they are the most forward ones -- true theta_lab 8-10 deg, where the (p,t) cross section is
# largest. That distorts the angular distribution, not just the normalisation.
# The default here is therefore wider. THE ACCEPTANCE SIMULATION MUST USE THE SAME VALUES, or data
# is corrected by an acceptance computed under a different cut.
THMIN="${THMIN:-5.0}"
THMAX="${THMAX:-175.0}"
GATE="${GATE:-$HERE/pid/triton_14C.json}"
# THE GATED FILE NEEDS THE PARAMETER-CONTAINER OBJECTS COPIED FROM A REFERENCE RECO FILE. Without
# them FairRunAna cannot initialise and the fit produces an empty output that ROOT opens without
# complaint -- "fRtdb->initContainers failed" in the log is the only sign. These reco files carry
# cbmout, BranchList and TimeBasedBranchList but no FileHeader, so cbmout is what to check for.
FREF="${FREF:-/mnt/h/a1954_C14_reco_hdb/}"
# FALLBACK REFERENCE. The primary directory does not hold every run -- run_0060 is absent from it.
# cbmout is a TFolder describing the branch structure, not per-run data, so any reco file of the
# same production serves. Without a fallback one missing reference silently drops a whole run from
# the analysis, which is the kind of loss that shows up later as an unexplained normalisation.
FREF2="${FREF2:-/mnt/f/a1954_C14_gf_xtr/}"
GEONAME="${GEONAME:-ATTPC_H300torr_RT}"
CACHETAG="${CACHETAG:-pt}"
LOG="$OUTDIR/logs"; mkdir -p "$GATEDIR" "$LOG"
source /home/yassid/fair_install/FairSoft/install/bin/thisroot.sh 2>/dev/null
set +u; source "$REPO/build/config.sh" >/dev/null 2>&1; set -u
export ROOT_INCLUDE_PATH="$REPO/build/include:$HOME/fair_install/FairRootInstall/include"
[ -f "$GATE" ] || { echo "no triton gate at $GATE"; exit 1; }
cd "$HERE"

one() {
  r=$1
  L="$LOG/${r}.log"; : > "$L"
  # 1. gate: IC 950-1350 + single pulse on the event, triton polygon + backward on the track
  ref="$FREF${r}_reco.root"
  [ -s "$ref" ] || ref="$FREF2${r}_reco.root"
  [ -s "$ref" ] || { echo "$r NO_REFERENCE (tried $FREF and $FREF2)"; return 1; }
  if [ ! -s "$GATEDIR/${r}_reco.root" ]; then
    root -b -q -l "pipeline/gate_events_C14.C(\"$r\",\"$IN\",\"$GATEDIR\",\"$ref\",950,1350,1050,1250,2.85,\"$GATE\",90.0,200,800,1500)" >> "$L" 2>&1
  fi
  [ -s "$GATEDIR/${r}_reco.root" ] || { echo "$r GATE_FAILED"; return 1; }
  # Fail HERE, not three stages later. Without the parameter-container objects copied from the
  # reference reco file, FairRunAna cannot initialise and the fit writes an EMPTY genfit file that
  # ROOT opens without complaint -- the only sign is "fRtdb->initContainers failed" in a log.
  # The key to check is cbmout: these reco files carry cbmout, BranchList and TimeBasedBranchList
  # but NO FileHeader, so testing for FileHeader rejects perfectly good input.
  root -b -q -l -e "TFile*f=TFile::Open(\"$GATEDIR/${r}_reco.root\"); \
    printf(\"HDRCHK %d\\n\", (f && !f->IsZombie() && f->Get(\"cbmout\")) ? 1 : 0);" 2>/dev/null \
    | grep -q "HDRCHK 1" || { echo "$r NO_cbmout -- check FREF=$FREF"; return 1; }
  # 2. fit: GENFIT + CATIMA, material effects ON, TRITON
  if [ ! -s "$OUTDIR/${r}_genfit.root" ]; then
    root -b -q -l "pipeline/fitGenfit_C14.C(\"$r\",-1,\"$GATEDIR\",\"\",\"$OUTDIR\",-2.85,2,5,\"\",4.0,$THMIN,$THMAX,kTRUE,kFALSE,\"triton\",\"$GEONAME\",kFALSE,kTRUE,0.0,1,kTRUE,kFALSE,kFALSE,kFALSE)" >> "$L" 2>&1
  fi
  [ -s "$OUTDIR/${r}_genfit.root" ] || { echo "$r FIT_FAILED"; return 1; }
  grep -q "dE/dx from CATIMA" "$L" || { echo "$r CATIMA_NOT_ENABLED"; return 1; }
  echo "$r done"
}
export -f one; export IN GATEDIR OUTDIR GATE GEONAME LOG HERE FREF FREF2 THMIN THMAX
printf "%s\n" $RUNS | xargs -P "$NPAR" -I{} bash -c 'one "$@"' _ {}

echo
echo "== kinematics cache: triton ejectile, 12C residual, theta window [$THMIN,$THMAX] =="
CSV=$(echo $RUNS | tr ' ' ',')
root -b -q -l "pp/ex_C14.C(\"$CSV\",\"$OUTDIR\",159.75,5.0,\"$CACHETAG\",3.016049,12.000000,\"(p,t)\",\"genfit\",kTRUE)" 2>&1 | tail -6
