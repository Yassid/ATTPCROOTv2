#!/usr/bin/env bash
# Browser explorer for 16C(p,d)15C (2026-11-14): the INTEGRAL-GATE selection beside the adopted one.
#
#     integral_gate   H2 float-tau re-reco (HDBSCAN + PSAMultiFit float tau) + genfit/CATIMA, deuteron
#                     gate on the trace-integral plane (gate_pd_h2_integral.json), 931 genfitter
#                     duplicate rows removed      /mnt/f/a1975/caches/pd_h2int/pd_kin_h2int_dedup.root
#     adopted         the cache behind the adopted (p,d) result (TC + PSAMax, peak-amplitude gate)
#                                                 /home/yassid/a1975_analysis/pd/data/pd_kin.root
# Both: same producer (cache_pd_run.C), same columns ke/theta, IC 950-1350. Finder, PSA AND gate all
# differ, so a difference here is not the gate alone. EBEAM defaults to 183.5 = run_pd_chain.sh.
# Structure copied from open_explorer_pd_catima.sh (see its notes on columns and collapsed fits).
#
#   ./open_explorer_pd_h2int.sh [Ebeam]
set -o pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EBEAM="${1:-183.5}"
TMP="${PDTMP:-/mnt/f/a1975/caches/.explorer}"
OUT="$HOME/a1975_C16_pd_integral_gate_explorer.html"
MK="$HERE/../../../a2091/UKF/pp/make_explorer_html.C"

C_CAT="/mnt/f/a1975/caches/pd_h2int/pd_kin_h2int_dedup.root"
C_OFF="/home/yassid/a1975_analysis/pd/data/pd_kin.root"
mkdir -p "$TMP"
set +u; source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh" >/dev/null 2>&1; set -u

for c in "$C_CAT" "$C_OFF"; do
  [[ -s "$c" ]] || { echo "ERROR: missing cache $c"; exit 1; }
done
[[ -s "$MK" ]] || { echo "ERROR: missing $MK"; exit 1; }

root -b -l -q "$HERE/mkexp_pp.C(\"$C_CAT\",\"$TMP/exp_pd_h2int.root\",1e8,950,1350,\"ke\",\"theta\")"
root -b -l -q "$HERE/mkexp_pp.C(\"$C_OFF\",\"$TMP/exp_pd_adopted.root\",1e8,950,1350,\"ke\",\"theta\")"

# 16C(p,d)15C: beam 16C, target p, ejectile d, residual 15C, beamA=16.
# 15C levels: g.s., the 0.740 5/2+ (the only bound excited state), Sn at 1.218, then unbound.
root -b -l -q "$MK(\"$TMP/exp_pd_h2int.root\",\"$OUT\",\"16C(p,d)15C  integral gate vs adopted\",$EBEAM,16.0147013,1.00782503,2.0135532,15.0105993,16,\"0:g.s.,0.740:5/2+,1.218:Sn,3.103,4.220,4.657\",\"$TMP/exp_pd_adopted.root\",\"\",\"integral_gate,adopted\")"
[[ -s "$OUT" ]] || { echo "ERROR: explorer not written"; exit 1; }
python3 "$HERE/add_keoff.py" "$OUT"

# Claude Code has no rendering surface, so the page has to be handed to a Windows browser.
if grep -qi microsoft /proc/version 2>/dev/null; then
  cp "$OUT" /mnt/c/Users/Yassid/Desktop/ 2>/dev/null && echo "copied to the Windows Desktop"
  WINPATH="$(wslpath -w "/mnt/c/Users/Yassid/Desktop/$(basename "$OUT")" 2>/dev/null)"
  [[ -n "$WINPATH" ]] && cmd.exe /c start "" "$WINPATH" >/dev/null 2>&1 && echo "opened in the default browser"
fi
echo "explorer -> $OUT"
