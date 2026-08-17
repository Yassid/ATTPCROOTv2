#!/usr/bin/env python3
"""Add a KE-offset control to a generated explorer page.

WHY: if GENFIT's energy (or angle) is biased, Ex is biased with it. A live offset on KE
lets you read off how large that bias would have to be to put the states on their lines --
without refitting anything. The offset is applied to the SAME ke everywhere: the KE cut,
the theta correction, the KE-vs-theta maps and the kinematics, so the page stays self
consistent and what you see cut is what you see plotted.

  python add_keoff.py <explorer.html>
"""
import sys, re
p = sys.argv[1]
s = open(p, encoding="utf-8").read()
orig = len(s)

# 1. register the control so readState() picks it up
s = s.replace("'kcDenom','kcPivot'];", "'kcDenom','kcPivot','keOff'];", 1)

# 2. the control itself, right after the beam group
anchor = '<button class="primary" id="zeroBtn">'
blk = ('''<div class="row"><label for="keOff">KE offset [MeV]</label><input type="number" id="keOff" step="0.05" value="0"></div>
      <input type="range" id="keOffR" min="-5" max="5" step="0.05" value="0">
      <div class="row"><span class="mono num" style="font-size:11.5px;color:var(--ink-3)">added to every KE: cut, maps and kinematics</span></div>
      ''')
s = s.replace(anchor, blk + anchor, 1)

# 3. apply it everywhere ke is read per track
s = s.replace("return d.c2[i] <= s.chi2\n      && thCorr(s, d.th[i], d.ke[i]) >= s.thLo && thCorr(s, d.th[i], d.ke[i]) <= s.thHi\n      && d.ke[i] >= s.keLo && d.ke[i] <= s.keHi",
              "const _k = d.ke[i] + (s.keOff||0);\n  return d.c2[i] <= s.chi2\n      && thCorr(s, d.th[i], _k) >= s.thLo && thCorr(s, d.th[i], _k) <= s.thHi\n      && _k >= s.keLo && _k <= s.keHi", 1)
n1 = s.count("const ke = d.ke[i], th = thCorr(s, d.th[i], ke);")
s = s.replace("const ke = d.ke[i], th = thCorr(s, d.th[i], ke);",
              "const ke = d.ke[i] + (s.keOff||0), th = thCorr(s, d.th[i], ke);")
# the comparison-overlay loop in exOnly() -- different shape, patched separately or the
# overlay would silently keep using the uncorrected KE while the main histogram moved
_o = "    const [ex] = kine2b(s.ebeam, thCorr(s, d.th[i], d.ke[i])*Math.PI/180, d.ke[i]);"
_n = ("    const _k = d.ke[i] + (s.keOff||0);\n"
      "    const [ex] = kine2b(s.ebeam, thCorr(s, d.th[i], _k)*Math.PI/180, _k);")
n2 = s.count(_o)
s = s.replace(_o, _n)

# 4. keep slider and number box in step, and redraw
hook = """
<script>
/* The page redraws through render(), wired to the 'input' event for every id in IDS.
   keOff is in IDS so the number box is already live; the slider is not, so it is wired
   here. Call render() directly and fall back to an 'input' event -- NOT 'change', which
   nothing listens for. */
(function(){
  const n=document.getElementById('keOff'), r=document.getElementById('keOffR');
  if(!n||!r) return;
  const fire=()=>{ if(typeof render==='function') render();
                   else n.dispatchEvent(new Event('input',{bubbles:true})); };
  r.addEventListener('input', ()=>{ n.value=r.value; fire(); });
  n.addEventListener('input', ()=>{ r.value=n.value; });
  n.addEventListener('change',()=>{ r.value=n.value; fire(); });
})();
</script>
"""
s = s.replace("</body>", hook + "</body>", 1)
s0 = s

# ---------------------------------------------------------------------------
# Make the states of interest VISIBLE on the KE-vs-theta panels.
# The stock palette gives distinct colours only to the first two reference levels;
# everything from index 2 up is grey (--refdim) and dashed, which on a viridis heat map
# is invisible. For 15C the level that matters most after the g.s./0.740 doublet is
# 3.103, and it was being drawn in grey dashes.
_old_col = "const refColor = i => css(i === 0 ? '--accent' : i === 1 ? '--ref2' : '--refdim');"
_new_col = ("const REF_COL = ['#ff7f0e','#2ca02c','#8f9aa4','#e6194b','#911eb4','#8f9aa4'];\n"
            "const refColor = i => REF_COL[i] || css('--refdim');")
if _old_col in s0:
    s0 = s0.replace(_old_col, _new_col, 1)
    # solid for g.s., 0.740 and 3.103; dashed for the Sn threshold and the higher ones
    s0 = s0.replace("const loci = REF_EX.map((r, i) => [r.ex, refColor(i), i < 2 ? [] : [5,4]]);",
                    "const loci = REF_EX.map((r, i) => [r.ex, refColor(i), (i === 2 || i > 3) ? [5,4] : []]);", 1)
    print("  reference loci recoloured: g.s. orange, 0.740 green, 3.103 red (solid)")
else:
    print("  WARNING: refColor not found, loci left as they were")

# ---------------------------------------------------------------------------
# Give the Ex-vs-theta_cm panel its own axis limits. As generated it has NO zoom at all:
# theta_cm is hardcoded 0..180 and Ex is tied to the main spectrum's exLo/exHi, so the only
# way to magnify a region was to rescale the Ex histogram too. cmLo/cmHi are slice markers,
# not limits. Four new controls, defaulting to the old behaviour.
_ids = "'kcDenom','kcPivot','keOff'];"
if _ids in s0:
    s0 = s0.replace(_ids, "'kcDenom','kcPivot','keOff','etCmLo','etCmHi','etExLo','etExHi'];", 1)
_anchor = '<div class="row"><label for="etExBins">'
_ctrl = ('<div class="row"><label>E<sub>x</sub> vs &theta;<sub>cm</sub> zoom</label></div>\n'
         '      <div class="pair"><input type="number" id="etCmLo" step="5" value="0" title="theta_cm min">'
         '<input type="number" id="etCmHi" step="5" value="180" title="theta_cm max"></div>\n'
         '      <div class="pair"><input type="number" id="etExLo" step="0.5" value="-5" title="Ex min">'
         '<input type="number" id="etExHi" step="0.5" value="10" title="Ex max"></div>\n      ')
if _anchor in s0:
    s0 = s0.replace(_anchor, _ctrl + _anchor, 1)
# bin on the panel's own ranges
_oldw = "const etW = (s.exHi - s.exLo)/s.etExBins;"
_neww = "const etW = ((s.etExHi ?? s.exHi) - (s.etExLo ?? s.exLo))/s.etExBins;"
s0 = s0.replace(_oldw, _neww, 1)
_oldf = "    const cx = Math.floor(tcm/180*s.cmBins), be = Math.floor((ex - s.exLo)/etW);"
_newf = ("    const _cl = s.etCmLo ?? 0, _ch = s.etCmHi ?? 180, _el = s.etExLo ?? s.exLo;\n"
         "    const cx = Math.floor((tcm - _cl)/(_ch - _cl)*s.cmBins), be = Math.floor((ex - _el)/etW);")
s0 = s0.replace(_oldf, _newf, 1)
# and draw on them
_oldd = "views.cEt = drawHeat($('cEt'), r.hEt, s.cmBins, s.etExBins, 0, 180, s.exLo, s.exHi,"
_newd = ("views.cEt = drawHeat($('cEt'), r.hEt, s.cmBins, s.etExBins, (s.etCmLo ?? 0), (s.etCmHi ?? 180),"
         " (s.etExLo ?? s.exLo), (s.etExHi ?? s.exHi),")
if _oldd in s0:
    s0 = s0.replace(_oldd, _newd, 1)
    print("  Ex vs theta_cm panel: independent zoom on both axes")
else:
    print("  WARNING: Et draw call not found")

s = s0
open(p, "w", encoding="utf-8").write(s)
print(f"  patched {p}: +{len(s)-orig} bytes, compute loops touched {n1}+{n2}")
