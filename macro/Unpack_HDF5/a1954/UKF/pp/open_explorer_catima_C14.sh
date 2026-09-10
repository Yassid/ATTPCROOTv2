#!/usr/bin/env bash
# Browser explorer for 14C(p,p') on the ADOPTED production: GENFIT + CATIMA material effects,
# chi2/ndf < 5 (plots/proton_kin_cat5.root, uncorrected theta). Self-contained, so it can be
# handed to anyone -- one HTML file, no ROOT needed.
#
# Built with the a2091 generator/template, not the a1954 one next to this script: that newer page
# has the live theta-correction block, whose formula is exactly apply_theta_corr_C14.C
#     theta' = theta - (360/N) * (KE - pivot)
# so the angle-corrected caches are reproduced in the page rather than baked in:
#     N = 2769.2308, pivot 0   ->  0.130 deg/MeV  = cat5_s013 (the document's adopted value)
#     N = 6428.5714, pivot 0   ->  0.056 deg/MeV  = cat5_tc   (the one that keeps the inelastic flat)
#
# The page opens on the adopted working point: Ebeam 159.75 (anchored on 6.094), correction ON at
# 0.130 deg/MeV, vertex slab 10-490 mm. All of them are live controls.
#
#   ./open_explorer_catima_C14.sh [Ebeam]
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
EBEAM="${1:-159.75}"
CACHE="$HERE/plots/proton_kin_cat5.root"
MK="$HERE/../../../a2091/UKF/pp/make_explorer_html.C"
OUT="$HOME/a1954_C14_pp_catima_explorer.html"
[[ -s "$CACHE" ]] || { echo "ERROR: missing cache $CACHE"; exit 1; }
[[ -s "$MK" ]] || { echo "ERROR: missing $MK"; exit 1; }

set +u; source "$HOME/fair_install/ATTPCROOTv2-OpenKF/build/config.sh" >/dev/null 2>&1; set -u

# 14C levels: g.s., the 6-7 MeV group, the isolated 8.317 2+, and Sn (above it is continuum)
LEVELS="0:g.s.,6.094:1-,6.728:3-,7.012:2+,7.341:2-,8.317:2+,8.176:Sn"
root -b -l -q "$MK(\"$CACHE\",\"$OUT\",\"14C(p,p')\",$EBEAM,14.003242,1.007825,1.007825,14.003242,14,\"$LEVELS\",\"\",\"\",\"GENFIT_CATIMA\")" || true
[[ -s "$OUT" ]] || { echo "ERROR: explorer not written"; exit 1; }

# the generator stamps a2091; and set the adopted defaults on the controls
python3 - "$OUT" <<'EOF'
import re, sys
p = sys.argv[1]
s = open(p, encoding='utf-8').read()
n0 = len(s)
subs = [
    ('"eyebrow":"a2091 . AT-TPC"', '"eyebrow":"a1954 . AT-TPC . GENFIT+CATIMA"'),
    ('id="hdrEyebrow">a2091 &middot; AT-TPC', 'id="hdrEyebrow">a1954 &middot; AT-TPC'),
    ('<input type="checkbox" id="kcOn">', '<input type="checkbox" id="kcOn" checked>'),
    ('id="kcDenom" step="50" value="4000"', 'id="kcDenom" step="50" value="2769.2308"'),
    ('id="kcSlopeR" min="-0.5" max="0.5" step="0.002" value="0.09"', 'id="kcSlopeR" min="-0.5" max="0.5" step="0.002" value="0.13"'),
    ('id="kcPivot" step="0.5" value="13"', 'id="kcPivot" step="0.5" value="0"'),
    ('id="kcPivotR" min="0" max="60" step="0.5" value="13"', 'id="kcPivotR" min="0" max="60" step="0.5" value="0"'),
    ('id="vzLo" step="10" value="-100"', 'id="vzLo" step="10" value="10"'),
    ('id="vzHi" step="10" value="1100"', 'id="vzHi" step="10" value="490"'),
    # the template presets the pivot to the data's median KE (2.5 MeV here) at load, overriding
    # the value above; the caches were made about pivot 0, so a 0.33 deg offset would creep in
    ("      $('kcPivot').value = (Math.round(med*2)/2).toFixed(1);\n"
     "      $('kcNote').textContent = 'pivot preset to median KE = ' + med.toFixed(1) + ' MeV';",
     "      $('kcNote').textContent = 'pivot 0, as in apply_theta_corr_C14 (median KE = ' + med.toFixed(1) + ' MeV)';"),
]
for a, b in subs:
    if s.count(a) != 1:
        sys.exit(f"ERROR: expected exactly one '{a}' in the page, found {s.count(a)}")
    s = s.replace(a, b)
open(p, 'w', encoding='utf-8').write(s)
# the zero-data trap: a page can be written with an empty set and still report success
m = re.search(r'"GENFIT_CATIMA":\{"ke":\[([^\]]*)\]', s)
nke = len(m.group(1).split(',')) if m and m.group(1) else 0
print(f"patched {len(subs)} defaults; GENFIT_CATIMA carries {nke} tracks")
if nke == 0:
    sys.exit("ERROR: the page carries no data")
EOF

cp "$OUT" /mnt/c/Users/Yassid/Desktop/ && echo "copied to the Windows Desktop"
echo "explorer -> $OUT"
