#!/usr/bin/env python3
"""Add the per-run Ex map to a generated C15d explorer page.

    python3 dd/inject_run_panel.py <explorer.html> [out.html]

WHY A MAP AND NOT A BAR CHART. The question this panel exists to answer is "does this structure
belong to one run or to all of them?" -- a per-run count cannot answer that, because a run can
contribute plenty of tracks and still put them at a different Ex. So the panel draws run (rows)
against Ex (columns), colour = counts, normalised PER ROW so a long run and a short one are
directly comparable. A feature present in every row is physics; a feature confined to a band of
rows is an artefact of those runs.

The panel calls the page's OWN kine2b/thCorr/pass at draw time rather than reimplementing the
two-body kinematics. The page's script is at global scope, so those are reachable; a second copy
would be a second thing to keep in step with the beam energy and the cuts.

The run cut itself lives in the page's cut model (runLo/runHi in IDS, applied in pass()), so
dragging here writes those inputs exactly as the IC panel writes icLo/icHi.
"""
import re
import sys


def find_grid_close(lines):
    start = next(i for i, l in enumerate(lines) if 'class="grid"' in l)
    depth = 0
    for i in range(start, len(lines)):
        depth += len(re.findall(r"<div\b", lines[i]))
        depth -= len(re.findall(r"</div>", lines[i]))
        if depth == 0:
            return i
    raise SystemExit("could not find the end of the .grid container")


CARD = '''
      <div class="card wide" id="cardRun">
        <h2>E<sub>x</sub> per run &mdash; drag vertically to select runs</h2>
        <div class="legend">
          <span>each row is one run, normalised to its own maximum</span>
          <span class="spacer" style="flex:1"></span>
          <button id="runAll" style="font-size:11px;padding:2px 8px">all runs</button>
        </div>
        <div class="canvas-wrap"><canvas id="cRun"></canvas><div class="tip" id="tipRun"></div></div>
        <div class="legend" id="runNote"></div>
      </div>
'''

SCRIPT = '''
<script>
/* ---- per-run Ex map (injected by dd/inject_run_panel.py) --------------------------------- */
(function(){
  const cv = document.getElementById('cRun'), tip = document.getElementById('tipRun');
  const note = document.getElementById('runNote');
  if (!cv || typeof SETS === 'undefined' || !HAS_RUN) return;
  const css = (n) => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
  const loI = document.getElementById('runLo'), hiI = document.getElementById('runHi');

  const d0 = SETS[ACTIVE];
  const RUNS = Array.from(new Set(Array.from(d0.run))).sort((a,b) => a-b);
  const IDX = new Map(RUNS.map((r,i) => [r,i]));

  let drag = null, geo = null;

  const setRuns = (a,b) => {
    if (!loI || !hiI) return;
    loI.value = a; hiI.value = b;
    loI.dispatchEvent(new Event('input', {bubbles:true}));
    hiI.dispatchEvent(new Event('input', {bubbles:true}));
    refresh();
  };

  /* Build the map with the RUN CUT LIFTED, so the picture always shows every run and the
     selection is drawn over it. A map that hid the runs you had deselected could not be used
     to decide which runs to select. */
  function build(){
    const s = readState();
    const sAll = Object.assign({}, s, {runLo: -1e9, runHi: 1e9});
    const nB = 90, lo = s.exLo, hi = s.exHi, w = (hi-lo)/nB;
    const d = SETS[ACTIVE];
    const grid = RUNS.map(() => new Float64Array(nB));
    const tot  = new Float64Array(RUNS.length);
    for (let i = 0; i < d.ke.length; ++i){
      if (!pass(d, i, sAll)) continue;
      const [ex] = kine2b(sAll.ebeam, thCorr(sAll, d.th[i], d.ke[i])*Math.PI/180, d.ke[i]);
      if (!isFinite(ex)) continue;
      const b = Math.floor((ex-lo)/w); if (b < 0 || b >= nB) continue;
      const r = IDX.get(d.run[i]); if (r === undefined) continue;
      grid[r][b]++; tot[r]++;
    }
    return {grid, tot, nB, lo, hi, w, s};
  }

  let M = null;
  function refresh(){ M = build(); draw(); }

  function draw(){
    if (!M) return;
    const dpr = window.devicePixelRatio || 1;
    const cw = cv.clientWidth || 900;
    const rowH = Math.max(3, Math.min(9, Math.floor(420/Math.max(1,RUNS.length))));
    const chh = RUNS.length*rowH + 46;
    cv.width = cw*dpr; cv.height = chh*dpr; cv.style.height = chh + 'px';
    const x = cv.getContext('2d'); x.setTransform(dpr,0,0,dpr,0,0);
    x.clearRect(0,0,cw,chh);
    const L = 46, R = 12, T = 8, B = 30, pw = cw-L-R, ph = RUNS.length*rowH;
    geo = {L,T,pw,rowH,ph};

    const cols = [css('--plot'), css('--data')];
    for (let r = 0; r < RUNS.length; ++r){
      const row = M.grid[r];
      let mx = 0; for (let b = 0; b < M.nB; ++b) if (row[b] > mx) mx = row[b];
      for (let b = 0; b < M.nB; ++b){
        const v = mx > 0 ? row[b]/mx : 0;
        if (v <= 0) continue;
        x.globalAlpha = 0.08 + 0.92*Math.sqrt(v);      // sqrt so a weak line still reads
        x.fillStyle = cols[1];
        x.fillRect(L + b*pw/M.nB, T + r*rowH, pw/M.nB + 0.6, rowH - 0.4);
      }
    }
    x.globalAlpha = 1;

    // selected run band
    const a = parseFloat(loI.value), b2 = parseFloat(hiI.value);
    const sel = drag ? [Math.min(drag.a,drag.b), Math.max(drag.a,drag.b)]
                     : (isFinite(a) && isFinite(b2) ? [a,b2] : null);
    if (sel){
      const i0 = RUNS.findIndex(r => r >= sel[0]);
      let i1 = -1; for (let i = RUNS.length-1; i >= 0; --i) if (RUNS[i] <= sel[1]) { i1 = i; break; }
      if (i0 >= 0 && i1 >= i0){
        x.strokeStyle = css('--accent'); x.lineWidth = 2;
        x.strokeRect(L-2, T + i0*rowH - 1, pw+4, (i1-i0+1)*rowH + 2);
      }
    }

    // axes
    x.strokeStyle = css('--line'); x.lineWidth = 1;
    x.beginPath(); x.moveTo(L,T); x.lineTo(L,T+ph); x.lineTo(L+pw,T+ph); x.stroke();
    x.fillStyle = css('--ink-3');
    x.font = '10.5px ui-monospace,SFMono-Regular,"DejaVu Sans Mono",Menlo,monospace';
    x.textAlign = 'center';
    const step = (M.hi-M.lo) <= 12 ? 2 : 5;
    for (let e = Math.ceil(M.lo/step)*step; e <= M.hi; e += step){
      const px = L + (e-M.lo)/(M.hi-M.lo)*pw;
      x.fillText(e, px, chh-14);
      x.globalAlpha = .25; x.strokeStyle = css('--grid');
      x.beginPath(); x.moveTo(px,T); x.lineTo(px,T+ph); x.stroke(); x.globalAlpha = 1;
    }
    x.fillText('E_x [MeV]', L+pw/2, chh-2);
    x.textAlign = 'right';
    for (let r = 0; r < RUNS.length; r += Math.max(1, Math.round(RUNS.length/18)))
      x.fillText(RUNS[r], L-5, T + r*rowH + rowH - 0.5);
  }

  const runAt = (e) => {
    const rect = cv.getBoundingClientRect();
    const i = Math.floor((e.clientY - rect.top - geo.T)/geo.rowH);
    return RUNS[Math.max(0, Math.min(RUNS.length-1, i))];
  };

  cv.addEventListener('mousedown', (e) => { if (!geo) return; const r = runAt(e); drag = {a:r,b:r}; e.preventDefault(); });
  cv.addEventListener('mousemove', (e) => {
    if (!geo || !M) return;
    if (drag){ drag.b = runAt(e); draw(); }
    const rect = cv.getBoundingClientRect();
    const ri = Math.floor((e.clientY - rect.top - geo.T)/geo.rowH);
    if (ri < 0 || ri >= RUNS.length){ tip.style.display = 'none'; return; }
    const ex = M.lo + ((e.clientX - rect.left) - geo.L)/geo.pw*(M.hi-M.lo);
    tip.style.display = 'block';
    tip.style.left = Math.min(e.clientX - rect.left + 12, rect.width-140) + 'px';
    tip.style.top = Math.max(4, ri*geo.rowH - 10) + 'px';
    tip.innerHTML = 'run ' + RUNS[ri] + '<br>' + M.tot[ri].toLocaleString('en-US')
                  + ' tracks<br>E_x ' + ex.toFixed(2);
  });
  addEventListener('mouseup', () => {
    if (!drag) return;
    const a = Math.min(drag.a, drag.b), b = Math.max(drag.a, drag.b);
    drag = null; setRuns(a, b);
  });
  cv.addEventListener('mouseleave', () => { tip.style.display = 'none'; });

  const allBtn = document.getElementById('runAll');
  if (allBtn) allBtn.addEventListener('click', () => setRuns(RUNS[0], RUNS[RUNS.length-1]));
  [loI,hiI].forEach(el => el && el.addEventListener('input', draw));

  // the map depends on ebeam and every cut, so rebuild whenever the page re-renders
  const _render = window.render;
  window.render = function(){ _render.apply(this, arguments); refresh(); };

  addEventListener('resize', draw);
  matchMedia('(prefers-color-scheme:dark)').addEventListener('change', draw);
  new MutationObserver(draw).observe(document.documentElement, {attributes:true, attributeFilter:['data-theme']});

  note.textContent = RUNS.length + ' runs, ' + RUNS[0] + '-' + RUNS[RUNS.length-1]
    + ' \\u00b7 each row normalised to its own maximum, so run length does not set brightness';
  refresh();
})();
</script>
'''


def main():
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    src = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else src
    lines = open(src, encoding="utf-8").read().split("\n")
    if 'id="cardRun"' in "\n".join(lines):
        raise SystemExit("panel already present in " + src)
    close = find_grid_close(lines)
    lines[close:close] = CARD.strip("\n").split("\n")
    lines.append(SCRIPT)
    open(out, "w", encoding="utf-8").write("\n".join(lines))
    print("injected per-run Ex map -> %s" % out)


if __name__ == "__main__":
    main()
