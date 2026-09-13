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

# The 1n phase space baked by phasespace_1n_C14.C, embedded so the global-fit panel can use the
# exact shape the ROOT macros fit with (the panel also generates its own, live, at the page's E_beam).
PS_ROOT="$HERE/plots/phasespace_1n_C14.root"
PS_JSON="$HERE/plots/.phasespace_1n_C14.json"
[[ -s "$PS_ROOT" ]] || { echo "ERROR: missing $PS_ROOT -- run phasespace_1n_C14.C"; exit 1; }
root -b -l -q -e 'TFile f("'"$PS_ROOT"'"); auto h = (TH1D*)f.Get("hPS"); if (!h) { printf("NOHPS\n"); return; }
  printf("GFHPS{\"lo\":%g,\"hi\":%g,\"v\":[", h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());
  for (int b = 1; b <= h->GetNbinsX(); ++b) printf("%.6g%s", h->GetBinContent(b), b < h->GetNbinsX() ? "," : "");
  printf("]}\n");' 2>/dev/null | grep '^GFHPS' | sed 's/^GFHPS//' > "$PS_JSON" || true
[[ -s "$PS_JSON" ]] || { echo "ERROR: could not read hPS from $PS_ROOT"; exit 1; }

# the generator stamps a2091; and set the adopted defaults on the controls
python3 - "$OUT" "$HERE/explorer_globalfit_C14.js" "$PS_JSON" <<'EOF'
import re, sys
p = sys.argv[1]
s = open(p, encoding='utf-8').read()
n0 = len(s)
gf = open(sys.argv[2], encoding='utf-8').read()
if gf.count('/*@@GF_HPS@@*/null') != 1:
    sys.exit("ERROR: GF_HPS placeholder not found exactly once in explorer_globalfit_C14.js")
gf = gf.replace('/*@@GF_HPS@@*/null', open(sys.argv[3], encoding='utf-8').read().strip())
subs = [
    # ---- global-fit panel (explorer_globalfit_C14.js): render hook, exports, and the code itself
    ("  $('foot').textContent = CFG.tag", "  globalFit(s);\n  $('foot').textContent = CFG.tag"),
    ("const panelIds = () => ['cEx','cKT','cEt','cSl','cVz']", "const panelIds = () => ['cEx','cKT','cEt','cSl','cVz','cGf']"),
    ("cSl:'Ex_slice', cVz:'Ex_vs_vertexz'};", "cSl:'Ex_slice', cVz:'Ex_vs_vertexz', cGf:'Ex_globalfit'};"),
    ("    if (v.kind === 'hist') {", "    if (v.kind === 'gfit') { lines.push(...v.csv); continue; }\n    if (v.kind === 'hist') {"),
    # the stitched figure assumed exactly one wide panel at an even index; lay panels out in order,
    # wide ones at full width and their own aspect ratio
    ("""    const w = cs[0].width, h = cs[0].height, rows = Math.ceil(cs.length/2);
    const out = document.createElement('canvas');
    out.width = w*2; out.height = h*rows;
    const g = out.getContext('2d');
    g.fillStyle = css('--plot'); g.fillRect(0,0,out.width,out.height);
    cs.forEach((cv, i) => {
      if (cv.id === 'cVz') g.drawImage(cv, 0, Math.floor(i/2)*h, w*2, h);
      else g.drawImage(cv, (i%2)*w, Math.floor(i/2)*h, w, h);
    });""",
     """    const w = cs[0].width, h = cs[0].height, place = [];
    let y = 0, col = 0;
    for (const cv of cs) {
      if (cv.id === 'cVz' || cv.id === 'cGf') {
        if (col) { y += h; col = 0; }
        const hh = Math.round(cv.height/cv.width*w*2);
        place.push([cv, 0, y, w*2, hh]); y += hh;
      } else {
        place.push([cv, col*w, y, w, h]);
        if (col) { y += h; col = 0; } else col = 1;
      }
    }
    if (col) y += h;
    const out = document.createElement('canvas');
    out.width = w*2; out.height = y;
    const g = out.getContext('2d');
    g.fillStyle = css('--plot'); g.fillRect(0,0,out.width,out.height);
    for (const [cv, x0, y0, ww, hh] of place) g.drawImage(cv, x0, y0, ww, hh);"""),
    ("\ninit();\n</script>", "\n" + gf + "\ninit();\n</script>"),
    # ---- adopted defaults
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
