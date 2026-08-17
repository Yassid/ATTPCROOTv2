#!/usr/bin/env bash
# Rebuild the GenFit fork + AtReconstruction after the CATIMA dE/dx change, IN THE RIGHT ORDER.
#
# THE ORDER IS NOT OPTIONAL. ATTPCROOT compiles against the INSTALLED GenFit headers in
# $GENFIT/include/, which are COPIES made at install time -- not the source tree in
# $GENFIT/trackReps/include/. Rebuilding ATTPCROOT first therefore fails with
#
#     error: 'class genfit::MaterialEffects' has no member named 'setUseCatimaEnergyLoss'
#
# which looks like a broken patch and is really just a stale header copy.
#
# DO NOT RUN THIS WHILE A FIT PRODUCTION IS RUNNING. The jobs hold libgenfit2.so and
# libAtReconstruction.so mapped; relinking underneath them risks crashing them or, worse,
# splitting one production across two physics models -- the same blended-model failure that
# matFallback = kFALSE exists to prevent. The guard below refuses to start if fits are running.
#
#   ./rebuild_catima_eloss.sh [nproc]
set -uo pipefail

NPROC="${1:-8}"
GENFIT=/home/yassid/fair_install/GenFit
ATTPC=/home/yassid/fair_install/ATTPCROOTv2-OpenKF

# Count only ROOT processes whose command line names one of the fit/cache macros. A plain
# `pgrep -f` on those names also matches ANY shell whose own command line happens to contain the
# pattern -- including the one that launched this script -- so it can refuse for no reason, and a
# guard that cries wolf is a guard people start bypassing. Matching on the executable name
# (root.exe, never bash) makes a self-match impossible.
running=0
for p in $(pgrep -x root.exe 2>/dev/null); do
  [ -r "/proc/$p/cmdline" ] || continue
  # NOT `tr ... | grep -q`: grep -q exits on first match and closes the pipe, tr dies of SIGPIPE,
  # and `set -o pipefail` then reports the whole pipeline as FAILED -- so a matching process would
  # be counted as absent and the build would proceed during a live production. Whether tr notices
  # depends on the pipe buffer, so it fails intermittently. Read into a variable, then match.
  cmdline=$(tr '\0' ' ' < "/proc/$p/cmdline" 2>/dev/null || true)
  case "$cmdline" in
  *fitGenfitter* | *fitUKF* | *cache_pd_run* | *cache_pp_run*)
    running=$((running + 1))
    ;;
  esac
done
if [ "$running" -gt 0 ]; then
  echo "REFUSING: $running fit/cache process(es) still running. Relinking under them would"
  echo "corrupt the production. Wait for them to finish, or kill them deliberately first."
  exit 1
fi

# DO NOT use `make install`. This tree was never deployed that way: CMAKE_INSTALL_PREFIX is
# /usr/local (empty, and it would need root), while $GENFIT/lib/libgenfit2.so.2.2.0 is a
# BYTE-IDENTICAL hand copy of build_catima/lib/..., and $GENFIT/include/*.h are flat hand copies
# of the per-module source headers. `make install` would write somewhere nothing loads and leave
# the library the production actually uses untouched -- a build that "succeeds" and changes nothing.
echo "=== 0/3  backing up the library the production currently loads ==="
cp -a "$GENFIT/lib/libgenfit2.so.2.2.0" "$GENFIT/lib/libgenfit2.so.2.2.0.bak" \
  && echo "    -> lib/libgenfit2.so.2.2.0.bak"

echo "=== 1/3  GenFit (build only) ==="
make -C "$GENFIT/build_catima" -j"$NPROC" || { echo "GENFIT BUILD FAILED"; exit 1; }

echo "=== 2/3  deploying lib + changed headers by hand, the way this tree is wired ==="
cp "$GENFIT/build_catima/lib/libgenfit2.so.2.2.0" "$GENFIT/lib/libgenfit2.so.2.2.0" \
  || { echo "LIB COPY FAILED"; exit 1; }
cp "$GENFIT/trackReps/include/MaterialEffects.h" "$GENFIT/include/MaterialEffects.h" \
  || { echo "HEADER COPY FAILED"; exit 1; }

# Prove BOTH halves of the deploy actually landed before spending time on ATTPCROOT.
if ! grep -q "setUseCatimaEnergyLoss" "$GENFIT/include/MaterialEffects.h"; then
  echo "$GENFIT/include/MaterialEffects.h has no setUseCatimaEnergyLoss -- the header copy did"
  echo "not land. ATTPCROOT compiles against THIS copy, not the source tree, so stopping here."
  exit 1
fi
# Same SIGPIPE trap as the guard above: `nm | grep -q` reports FAILURE under pipefail precisely
# when grep succeeds early. Capture first, then match -- no pipe, no false verdict.
gfsyms=$(nm -D --defined-only "$GENFIT/lib/libgenfit2.so" 2>/dev/null || true)
case "$gfsyms" in
*catimaDEdx*) ;;
*)
  echo "$GENFIT/lib/libgenfit2.so has no catimaDEdx symbol -- the library copy did not land,"
  echo "or the build silently produced the non-CATIMA variant. Stopping before ATTPCROOT."
  exit 1
  ;;
esac
echo "    header and library both carry the new dE/dx path"

echo "=== 3/3  AtReconstruction ==="
cmake --build "$ATTPC/build" --target AtReconstruction -j"$NPROC" \
  || { echo "ATTPCROOT BUILD FAILED"; exit 1; }

echo
echo "=== built. Before trusting any physics from it: ==="
echo "  * the A/B must run on a SPREAD of runs, not one. The first CATIMA recommendation looked"
echo "    good because it was validated on run_0031 alone, which turned out to be 1 of only 6"
echo "    runs where matFX worked; the other 41 collapsed 40-94%."
echo "  * compare against the TABLE production with everything else identical, and use MEDIANS."
echo "  * SetCatimaELoss(kTRUE) replaces only the tabulated branch (beta*gamma < 0.05);"
echo "    SetCatimaELoss(kTRUE, kTRUE) also replaces Bethe-Bloch. Validate them separately."
