#!/usr/bin/env bash
# Browser explorer for 16C(p,d)15C (2026-11-15) with a LIVE IC window and the Ebeam(vz) knob.
#
#     tc_tau        H2 re-reco, TC finder + PSAMultiFit float tau, genfit/CATIMA, integral gate,
#                   686 genfitter duplicates removed   /mnt/f/a1975/caches/pd_tcint/pd_kin_tcint_dedup.root
#     adopted       the cache behind the adopted (p,d) result (TC + PSAMax, peak-amplitude gate)
#                                                      /home/yassid/a1975_analysis/pd/data/pd_kin.root
#     hdb_tau       HDBSCAN + float tau + integral gate (the 2026-11-14 test, ~11 % low)
#                                                      /mnt/f/a1975/caches/pd_h2int/pd_kin_h2int_dedup.root
#
# Unlike open_explorer_pd_h2int.sh this uses the a1954_Be12 generator: its page has run and IC
# controls, so the slim ntuples are written UNGATED with ic/run (mkexp_pp keepIcRun) and the IC
# window is set on the page. It opens at 950-1350 = run_pd_chain.sh. Collapsed fits (chi2 1e9)
# are still dropped by mkexp_pp.
# Ebeam(vz) knob: add_ebvz.py, default OFF, slope 10.9 MeV/m = CATIMA 16C 185 MeV in H2 300 torr.
#
#   ./open_explorer_pd_tctau.sh [Ebeam]
set -o pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EBEAM="${1:-183.5}"
TMP="${PDTMP:-/mnt/f/a1975/caches/.explorer}"
OUT="$HOME/a1975_C16_pd_tctau_explorer.html"
MK="$HERE/../../../a1954_Be12/UKF/pp/make_explorer_html.C"

C1="/mnt/f/a1975/caches/pd_tcint/pd_kin_tcint_dedup.root"
C2="/home/yassid/a1975_analysis/pd/data/pd_kin.root"
C3="/mnt/f/a1975/caches/pd_h2int/pd_kin_h2int_dedup.root"
mkdir -p "$TMP"
set +u; source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh" >/dev/null 2>&1; set -u

for c in "$C1" "$C2" "$C3" "$MK"; do
  [[ -s "$c" ]] || { echo "ERROR: missing $c"; exit 1; }
done

root -b -l -q "$HERE/mkexp_pp.C(\"$C1\",\"$TMP/expic_pd_tctau.root\",1e8,-1,1e9,\"ke\",\"theta\",true)"
root -b -l -q "$HERE/mkexp_pp.C(\"$C2\",\"$TMP/expic_pd_adopted.root\",1e8,-1,1e9,\"ke\",\"theta\",true)"
root -b -l -q "$HERE/mkexp_pp.C(\"$C3\",\"$TMP/expic_pd_h2int.root\",1e8,-1,1e9,\"ke\",\"theta\",true)"

# 16C(p,d)15C: 15C levels g.s., 0.740 5/2+ (only bound excited state), Sn 1.218, then unbound.
root -b -l -q "$MK(\"$TMP/expic_pd_tctau.root\",\"$OUT\",\"16C(p,d)15C  TC-tau vs adopted vs HDBSCAN-tau\",$EBEAM,16.0147013,1.00782503,2.0135532,15.0105993,16,\"0:g.s.,0.740:5/2+,1.218:Sn,3.103,4.220,4.657\",\"$TMP/expic_pd_adopted.root\",\"$TMP/expic_pd_h2int.root\",\"tc_tau,adopted,hdb_tau\")"
[[ -s "$OUT" ]] || { echo "ERROR: explorer not written"; exit 1; }

# open on the chain's beam window instead of the template's "everything" (-1..4000)
sed -i 's|id="icLo" step="10" value="-1"|id="icLo" step="10" value="950"|; s|id="icHi" step="10" value="4000"|id="icHi" step="10" value="1350"|' "$OUT"
grep -q 'id="icLo" step="10" value="950"' "$OUT" || echo "WARNING: IC default not set -- page opens ungated" >&2

python3 "$HERE/add_keoff.py" "$OUT"
python3 "$HERE/add_ebvz.py" "$OUT" 10.9 "CATIMA 16C 185 MeV in H2 300 torr: 10.9 MeV/m"
python3 "$HERE/add_exshift.py" "$OUT" -0.6309

# ADOPTED (p,d) DEFAULTS = a1975_analysis/pd/run_pd_chain.sh + data/explorer_params_16C_p_d_15C_
# adopted_cache_2026-08-20-12-15-11.json: theta correction ON, kcDenom 4864.9, kcPivot 15.5,
# keOff -0.1, EXSHIFT -0.6309 (exShift control above), chi2 < 5, vz 0-500, IC 950-1350, Ebeam 183.5
# (passed to the generator). The Reset button returns to these.
sed -i -e 's|id="kcDenom" step="50" value="[^"]*"|id="kcDenom" step="50" value="4864.9"|' \
       -e 's|id="kcPivot" step="0.5" value="[^"]*"|id="kcPivot" step="0.5" value="15.5"|' \
       -e 's|id="keOff" step="0.05" value="[^"]*"|id="keOff" step="0.05" value="-0.1"|' \
       -e 's|id="chi2" step="0.5" value="[^"]*"|id="chi2" step="0.5" value="5"|' \
       -e 's|id="vzLo" step="10" value="[^"]*"|id="vzLo" step="10" value="0"|' \
       -e 's|id="vzHi" step="10" value="[^"]*"|id="vzHi" step="10" value="500"|' \
       -e 's|<input type="checkbox" id="kcOn">|<input type="checkbox" id="kcOn" checked>|' "$OUT"
for want in 'id="kcDenom" step="50" value="4864.9"' 'id="kcPivot" step="0.5" value="15.5"' \
            'id="keOff" step="0.05" value="-0.1"' 'id="vzLo" step="10" value="0"' 'id="vzHi" step="10" value="500"' \
            'id="kcOn" checked' 'id="exShift" step="0.01" value="-0.6309"'; do
  grep -q "$want" "$OUT" || echo "WARNING: default not applied: $want" >&2
done

if grep -qi microsoft /proc/version 2>/dev/null; then
  cp "$OUT" /mnt/c/Users/Yassid/Desktop/ 2>/dev/null && echo "copied to the Windows Desktop"
  WINPATH="$(wslpath -w "/mnt/c/Users/Yassid/Desktop/$(basename "$OUT")" 2>/dev/null)"
  [[ -n "$WINPATH" ]] && cmd.exe /c start "" "$WINPATH" >/dev/null 2>&1 && echo "opened in the default browser"
fi
echo "explorer -> $OUT"
