#!/bin/bash
# Generate + digitise alpha+alpha events in parallel streams, for the trigger-efficiency
# statistics (see doc/Dec2014_alphas_manual.md section 5.5).
#
#   usage: bulk_sim.sh <nStreams> <eventsPerStream> [outRoot]
#
# Each stream runs in its OWN directory because both macros write ./data/attpcsim_in.root
# and ./data/attpcpar_in.root by name -- streams sharing a cwd would overwrite each other's
# files mid-run. Each also gets an explicit non-zero seed: the macro's default seed comes
# from time(NULL), which has one-second granularity, so streams launched together would
# otherwise produce byte-identical events.
#
# Reminder: AtTPC2Body alternates beam-only and reaction events, so N events give N/2
# reactions.
# NO `set -u` here: setup_fr19port.sh (and FairSoft's thisroot.sh underneath it) reference
# unset variables as a matter of course, so `set -u` kills the shell at the source line --
# silently, before anything is echoed.
NSTREAM=${1:-8}
NEV=${2:-5000}
OUT=${3:-/home/yassid/dec2014_sim_bulk}

source /home/yassid/fair_install/ATTPCROOTv2_fr19port/setup_fr19port.sh > /dev/null 2>&1
MACRODIR=$VMCWORKDIR/macro/Simulation/ATTPC/Dec2014_alphas
mkdir -p "$OUT"

echo "[bulk] $NSTREAM streams x $NEV events = $((NSTREAM*NEV)) events (~$((NSTREAM*NEV/2)) reactions)"
echo "[bulk] start $(date +%H:%M:%S)"

for i in $(seq 1 "$NSTREAM"); do
(
   d="$OUT/stream_$i"
   rm -rf "$d"; mkdir -p "$d/data"
   # rootlogon.C is read from the CURRENT directory and must load the dictionaries before
   # the macro is parsed -- without it every class comes out as "unknown type name".
   ln -sf "$MACRODIR/rootlogon.C" "$d/rootlogon.C"
   cd "$d" || exit 1

   seed=$((1000 + i * 7919))          # distinct, deterministic per stream
   timeout 7200 root -l -b -q "$MACRODIR/He4He4_sim_el.C($NEV,\"TGeant4\",$seed)" \
        > sim.log 2>&1
   rc1=$?
   if [ ! -s data/attpcsim_in.root ]; then
      echo "[s$i ] SIM FAILED rc=$rc1 -- see $d/sim.log"; exit 1
   fi

   # keepElectrons=kFALSE: the drifted-electron collection is ~90% of the volume and is
   # only needed for the truth-residual study, not for charge spectra.
   timeout 14400 root -l -b -q \
        "$MACRODIR/rundigi_sim.C(\"./data/attpcsim_in.root\",\"./data/digi.root\",kFALSE)" \
        > digi.log 2>&1
   rc2=$?
   if [ ! -s data/digi.root ]; then
      echo "[s$i ] DIGI FAILED rc=$rc2 -- see $d/digi.log"; exit 1
   fi

   # Per-event summed GetQHit(), the observable the trigger study consumes.
   timeout 3600 root -l -b -q \
        "$MACRODIR/../../../Unpack_GETDecoder2/dump_hits.C(\"./data/digi.root\",\"./qtot.txt\",\"qtot\",-1)" \
        > dump.log 2>&1
   n=$(wc -l < qtot.txt 2>/dev/null || echo 0)
   echo "[s$i ] done  events=$n  $(du -h data/digi.root 2>/dev/null | cut -f1)  $(date +%H:%M:%S)"
) &
done
wait

cat "$OUT"/stream_*/qtot.txt > "$OUT/qtot_sim_bulk.txt" 2>/dev/null
echo "[bulk] TOTAL events dumped: $(wc -l < "$OUT/qtot_sim_bulk.txt")"
echo "[bulk] end $(date +%H:%M:%S)"
echo "BULK_DONE"
