#!/usr/bin/env python3
"""Add the ion-chamber panel to a generated C15p explorer page.

    python3 pp/inject_ic_panel.py <explorer.html> <ic.json> [out.html]

The explorer template is shared with the other channels, so the IC panel is INJECTED into the
generated page rather than added to the template -- that keeps a C15p-only panel out of every
other analysis's viewer. The panel is appended as the last card in the existing .grid, styled
purely through the page's own CSS custom properties, so it follows the light/dark toggle with
no styling of its own to keep in sync.

ic.json carries {lo, hi, nb, gateLo, gateHi, nTot, nNoIC, nSP, all[], sp[]} as written by the
IC dump step: `all` is every track with a usable IC value, `sp` the single-pulse subset, which
is the one the gate is defined on (pid/ic_C15p.json records singlePulse: true).
"""
import json
import re
import sys


def find_grid_close(lines):
    """Index of the </div> that closes <div class="grid">, by depth counting."""
    start = next(i for i, l in enumerate(lines) if 'class="grid"' in l)
    depth = 0
    for i in range(start, len(lines)):
        depth += len(re.findall(r"<div\b", lines[i]))
        depth -= len(re.findall(r"</div>", lines[i]))
        if depth == 0:
            return i
    raise SystemExit("could not find the end of the .grid container")


CARD = '''
      <div class="card wide" id="cardIC">
        <h2>ion chamber &mdash; drag to select the beam window</h2>
        <div class="legend">
          <span><i class="swatch" style="background:var(--refdim)"></i>all multiplicities</span>
          <span><i class="swatch" style="background:var(--data)"></i>selected multiplicity</span>
          <span><i class="swatch" style="background:var(--accent);opacity:.35;height:9px;width:9px;border-radius:2px"></i>selection</span>
          <span class="spacer" style="flex:1"></span>
          <button id="icSnap" style="font-size:11px;padding:2px 8px">use [931, 1413]</button>
          <button id="icOff" style="font-size:11px;padding:2px 8px">clear</button>
          <label class="check" style="gap:5px"><input type="checkbox" id="icLog" checked> log y</label>
        </div>
        <div class="canvas-wrap"><canvas id="cIC" style="cursor:ew-resize"></canvas><div class="tip" id="tipIC"></div></div>
        <div class="legend" id="icNote"></div>
      </div>
'''

SCRIPT = '''
<script>
/* ---- ion-chamber panel (injected by pp/inject_ic_panel.py) -------------------------------
   Draws its own canvas and reads colour from the page's CSS custom properties at paint time,
   so the theme toggle and the OS preference both carry it without a second palette. */
(function(){
  const IC = /*__ICDATA__*/ null;
  if (!IC) return;
  const cv = document.getElementById('cIC'), tip = document.getElementById('tipIC');
  const logBox = document.getElementById('icLog'), note = document.getElementById('icNote');
  if (!cv) return;
  const W = (i) => IC.lo + (i + 0.5) * (IC.hi - IC.lo) / IC.nb;   // bin centre
  const css = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
  const fmt = (n) => n.toLocaleString('en-US');

  /* The panel does not own the selection -- the page's icLo/icHi inputs do, so the window the
     panel shades is by construction the window every other panel is filtered on. Dragging here
     writes those inputs and fires their 'input' event, which is what the page already listens to
     for a re-render. icLo < 0 means no window, matching apply_gate_C15p.C's convention. */
  const loI = document.getElementById('icLo'), hiI = document.getElementById('icHi');
  const win = () => {
    if (!loI || !hiI) return null;
    const lo = parseFloat(loI.value), hi = parseFloat(hiI.value);
    return (isFinite(lo) && lo >= 0 && isFinite(hi) && hi > lo) ? [lo, hi] : null;
  };
  const setWin = (lo, hi) => {
    if (!loI || !hiI) return;
    loI.value = Math.round(lo); hiI.value = Math.round(hi);
    loI.dispatchEvent(new Event('input', {bubbles: true}));
    hiI.dispatchEvent(new Event('input', {bubbles: true}));
    refresh();
  };

  const npI = document.getElementById('npLo'), npH = document.getElementById('npHi');
  /* byNp is indexed 0..3 for multiplicity 1,2,3,4+. The panel shades whatever the page's
     npLo/npHi currently select, so "what am I selecting" is answered on the plot itself. */
  const npSel = () => {
    const a = npI ? parseFloat(npI.value) : 1, b = npH ? parseFloat(npH.value) : 1;
    const lo = isFinite(a) ? a : 1, hi = isFinite(b) ? b : 1;
    const out = [];
    for (let m = 0; m < 4; ++m){ const mult = m + 1;           // index 3 means "4 or more"
      const inRange = (m < 3) ? (mult >= lo && mult <= hi) : (hi >= 4);
      if (inRange) out.push(m); }
    return out;
  };
  const sumOver = (idx) => {
    const h = new Float64Array(IC.nb);
    idx.forEach(m => { for (let i = 0; i < IC.nb; ++i) h[i] += IC.byNp[m][i]; });
    return h;
  };
  const ALL = sumOver([0,1,2,3]);

  let drag = null;   // {a, b} in IC units, while the mouse is down

  function refresh(){
    const w = drag ? [Math.min(drag.a, drag.b), Math.max(drag.a, drag.b)] : win();
    const idx = npSel(), SEL = sumOver(idx);
    let selTot = 0, inWin = 0;
    for (let i = 0; i < IC.nb; ++i){ selTot += SEL[i];
      if (w && W(i) >= w[0] && W(i) <= w[1]) inWin += SEL[i]; }
    const lbl = ['1','2','3','4+'].filter((_, m) => idx.includes(m)).join(',') || 'none';
    note.textContent = fmt(IC.nTot) + ' tracks \\u00b7 ' + fmt(IC.nNoIC) + ' with no usable IC \\u00b7 '
      + 'multiplicity ' + lbl + ' -> ' + fmt(selTot) + ' \\u00b7 '
      + (w ? ('window [' + Math.round(w[0]) + ', ' + Math.round(w[1]) + '] keeps ' + fmt(inWin)
              + ' = ' + (100 * inWin / Math.max(1, selTot)).toFixed(1) + '%')
           : 'no window \\u2014 drag across the plot to set one')
      + ' \\u00b7 all: ' + IC.cnt.map((c, m) => ['1','2','3','4+'][m] + ':' + fmt(c)).join('  ');
    draw();
  }

  let geo = null;
  function draw(){
    const dpr = window.devicePixelRatio || 1;
    const cw = cv.clientWidth || 900, chh = 260;
    cv.width = cw * dpr; cv.height = chh * dpr; cv.style.height = chh + 'px';
    const x = cv.getContext('2d'); x.setTransform(dpr, 0, 0, dpr, 0, 0);
    x.clearRect(0, 0, cw, chh);
    const L = 54, R = 12, T = 10, B = 30;
    const pw = cw - L - R, ph = chh - T - B;
    const lg = logBox.checked;
    const top = Math.max(1, ...ALL);
    const yv = (v) => lg ? (v <= 0 ? 0 : Math.log10(v + 1) / Math.log10(top + 1)) : v / top;
    const px = (w) => L + (w - IC.lo) / (IC.hi - IC.lo) * pw;
    const py = (v) => T + ph - yv(v) * ph;
    geo = {L, T, pw, ph, chh};

    // selection band first, so the histogram draws over it
    const w = drag ? [Math.min(drag.a, drag.b), Math.max(drag.a, drag.b)] : win();
    if (w){
      x.fillStyle = css('--accent-soft');
      x.fillRect(px(w[0]), T, px(w[1]) - px(w[0]), ph);
      x.strokeStyle = css('--accent'); x.lineWidth = 1; x.globalAlpha = .85;
      w.forEach(g => { x.beginPath(); x.moveTo(px(g), T); x.lineTo(px(g), T + ph); x.stroke(); });
      x.globalAlpha = 1;
    }
    // the window recorded in pid/ic_C15p.json, marked but never enforced, so a hand-dragged
    // selection can be compared against the one the batch macros use
    x.save();
    x.strokeStyle = css('--ref2'); x.lineWidth = 1; x.setLineDash([4, 3]); x.globalAlpha = .9;
    [IC.gateLo, IC.gateHi].forEach(g => {
      x.beginPath(); x.moveTo(px(g), T); x.lineTo(px(g), T + ph); x.stroke();
    });
    x.restore();

    // grid + axes
    x.strokeStyle = css('--grid'); x.lineWidth = 1;
    for (let g = 0; g <= 3000; g += 500){
      x.beginPath(); x.moveTo(px(g), T); x.lineTo(px(g), T + ph); x.stroke();
    }
    x.strokeStyle = css('--line'); x.beginPath();
    x.moveTo(L, T); x.lineTo(L, T + ph); x.lineTo(L + pw, T + ph); x.stroke();

    const bars = (arr, col, alpha) => {
      x.fillStyle = col; x.globalAlpha = alpha;
      const bw = pw / IC.nb;
      for (let i = 0; i < IC.nb; i++){
        if (arr[i] <= 0) continue;
        const h = T + ph - py(arr[i]);
        x.fillRect(px(W(i)) - bw / 2, py(arr[i]), Math.max(bw, .8), h);
      }
      x.globalAlpha = 1;
    };
    bars(ALL,              css('--refdim'), .40);
    bars(sumOver(npSel()), css('--data'),   .95);

    // labels
    x.fillStyle = css('--ink-3');
    x.font = '11px ui-monospace,SFMono-Regular,"DejaVu Sans Mono",Menlo,monospace';
    x.textAlign = 'center';
    for (let g = 0; g <= 3000; g += 500) x.fillText(g, px(g), chh - 10);
    x.fillText('IC max amplitude [ADC]', L + pw / 2, chh - 1);
    x.textAlign = 'right';
    x.fillText(lg ? 'log' : fmt(top), L - 6, T + 10);
    x.fillText('0', L - 6, T + ph);
    // the two cocktail components
    x.fillStyle = css('--ink-2'); x.textAlign = 'center';
    x.fillText('1168', px(1168), T + 11);
    x.fillText('2058', px(2058), T + 11);
  }

  const atX = (e) => {
    const r = cv.getBoundingClientRect();
    return IC.lo + ((e.clientX - r.left) - geo.L) / geo.pw * (IC.hi - IC.lo);
  };
  const clamp = (v) => Math.max(IC.lo, Math.min(IC.hi, v));

  cv.addEventListener('mousedown', (e) => {
    if (!geo) return;
    const v = clamp(atX(e));
    drag = {a: v, b: v};
    e.preventDefault();
  });
  cv.addEventListener('mousemove', (e) => {
    if (!geo) return;
    const r = cv.getBoundingClientRect(), mx = e.clientX - r.left;
    if (drag){ drag.b = clamp(atX(e)); refresh(); }
    const wv = atX(e);
    if (wv < IC.lo || wv > IC.hi){ tip.style.display = 'none'; return; }
    const i = Math.min(IC.nb - 1, Math.max(0, Math.floor((wv - IC.lo) / (IC.hi - IC.lo) * IC.nb)));
    tip.style.display = 'block';
    tip.style.left = Math.min(mx + 12, r.width - 150) + 'px';
    tip.style.top = '14px';
    tip.innerHTML = 'IC ' + Math.round(W(i)) + '<br>all ' + fmt(ALL[i])
                  + '<br>np 1 ' + fmt(IC.byNp[0][i]) + '<br>np 2 ' + fmt(IC.byNp[1][i]);
  });
  // A click without movement is not a zero-width window -- that would empty every panel and read
  // as a broken page. Anything under 15 ADC wide is treated as a click and clears the selection.
  addEventListener('mouseup', () => {
    if (!drag) return;
    const a = Math.min(drag.a, drag.b), b = Math.max(drag.a, drag.b);
    drag = null;
    if (b - a < 15){ if (loI) setWin(-1, 4000); else refresh(); return; }
    setWin(a, b);
  });
  cv.addEventListener('mouseleave', () => { tip.style.display = 'none'; });

  const snap = document.getElementById('icSnap'), off = document.getElementById('icOff');
  if (snap) snap.addEventListener('click', () => setWin(IC.gateLo, IC.gateHi));
  if (off)  off.addEventListener('click', () => { if (loI){ loI.value = -1; hiI.value = 4000;
              loI.dispatchEvent(new Event('input', {bubbles:true})); } refresh(); });
  [loI, hiI, npI, npH].forEach(el => el && el.addEventListener('input', refresh));

  logBox.addEventListener('change', draw);
  addEventListener('resize', draw);
  matchMedia('(prefers-color-scheme:dark)').addEventListener('change', draw);
  new MutationObserver(draw).observe(document.documentElement, {attributes:true, attributeFilter:['data-theme']});
  refresh();
})();
</script>
'''


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    src, icf = sys.argv[1], sys.argv[2]
    out = sys.argv[3] if len(sys.argv) > 3 else src
    lines = open(src, encoding="utf-8").read().split("\n")
    ic = json.load(open(icf, encoding="utf-8"))

    if 'id="cardIC"' in "\n".join(lines):
        raise SystemExit("panel already present in " + src)

    close = find_grid_close(lines)
    lines[close:close] = CARD.strip("\n").split("\n")

    script = SCRIPT.replace("/*__ICDATA__*/ null", json.dumps(ic, separators=(",", ":")))
    # after the page's own script, so nothing races the initial render
    lines.append(script)

    open(out, "w", encoding="utf-8").write("\n".join(lines))
    print("injected IC panel -> %s (%d bins, gate [%.0f, %.0f])"
          % (out, ic["nb"], ic["gateLo"], ic["gateHi"]))


if __name__ == "__main__":
    main()
