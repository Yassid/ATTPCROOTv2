#!/usr/bin/env python3
"""Build the ion-chamber trace viewer from a dump written by pid/dump_ic_traces_C15d.C.

    python3 pid/make_ic_viewer_C15d.py <ic_traces.json> <out.html>

The page shows one raw IC trace at a time with BOTH reductions drawn on it: the amplitude window
icsum_C15d.C actually uses ([1050,1250], shaded) and the wider window it counts pulses over
([800,1500], outlined), together with the pulse threshold and the two maxima those windows find.
The point is to make a disagreement between them visible on the pulse rather than inferable from
a histogram.

Palette is lifted verbatim from dd/explorer_template.html so the two pages read as one tool, and
carries the same light/dark handling: prefers-color-scheme for the OS, data-theme for the
viewer's own toggle.
"""
import json
import sys

TOKENS = """
:root{
  --bg:#f2f4f6; --panel:#ffffff; --panel-2:#f8f9fb;
  --line:#dadfe5; --line-2:#eaeef2; --grid:#e7ebef;
  --ink:#161b20; --ink-2:#4c5761; --ink-3:#7c8792;
  --accent:#a95c29; --accent-ink:#8d4c21; --accent-soft:#f4e6da;
  --data:#1d5f78; --fit:#c03f2b; --ref2:#2c7a58; --refdim:#8f9aa4;
  --plot:#ffffff; --good:#2c7a58; --warn:#a9761f;
  --shadow:0 1px 2px rgba(20,26,32,.05),0 6px 18px rgba(20,26,32,.05);
}
@media (prefers-color-scheme:dark){
  :root{
    --bg:#0d1115; --panel:#161b21; --panel-2:#1a2027;
    --line:#2a323a; --line-2:#222932; --grid:#232b33;
    --ink:#e7ebef; --ink-2:#a6b1bb; --ink-3:#78838d;
    --accent:#d9843f; --accent-ink:#e79a5c; --accent-soft:#2e2318;
    --data:#68b8d6; --fit:#ef6a52; --ref2:#4fb98a; --refdim:#7d8892;
    --plot:#10151a; --good:#4fb98a; --warn:#d9a441;
    --shadow:0 1px 2px rgba(0,0,0,.4),0 6px 18px rgba(0,0,0,.3);
  }
}
:root[data-theme="light"]{
  --bg:#f2f4f6; --panel:#ffffff; --panel-2:#f8f9fb;
  --line:#dadfe5; --line-2:#eaeef2; --grid:#e7ebef;
  --ink:#161b20; --ink-2:#4c5761; --ink-3:#7c8792;
  --accent:#a95c29; --accent-ink:#8d4c21; --accent-soft:#f4e6da;
  --data:#1d5f78; --fit:#c03f2b; --ref2:#2c7a58; --refdim:#8f9aa4;
  --plot:#ffffff; --good:#2c7a58; --warn:#a9761f;
  --shadow:0 1px 2px rgba(20,26,32,.05),0 6px 18px rgba(20,26,32,.05);
}
:root[data-theme="dark"]{
  --bg:#0d1115; --panel:#161b21; --panel-2:#1a2027;
  --line:#2a323a; --line-2:#222932; --grid:#232b33;
  --ink:#e7ebef; --ink-2:#a6b1bb; --ink-3:#78838d;
  --accent:#d9843f; --accent-ink:#e79a5c; --accent-soft:#2e2318;
  --data:#68b8d6; --fit:#ef6a52; --ref2:#4fb98a; --refdim:#7d8892;
  --plot:#10151a; --good:#4fb98a; --warn:#d9a441;
  --shadow:0 1px 2px rgba(0,0,0,.4),0 6px 18px rgba(0,0,0,.3);
}
"""

PAGE = """<title>ion-chamber traces &mdash; C15d</title>
<style>
__TOKENS__
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--ink);
     font:14px/1.5 ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,Helvetica,Arial,sans-serif}
.mono{font-family:ui-monospace,SFMono-Regular,"JetBrains Mono","DejaVu Sans Mono",Menlo,monospace}
header{display:flex;align-items:baseline;gap:14px;flex-wrap:wrap;
       padding:16px 22px 10px;border-bottom:1px solid var(--line)}
h1{margin:0;font-size:17px;font-weight:650;letter-spacing:-.01em}
.eyebrow{font-size:10.5px;letter-spacing:.09em;text-transform:uppercase;color:var(--ink-3)}
.spacer{flex:1}
button{background:var(--panel);color:var(--ink);border:1px solid var(--line);border-radius:7px;
       padding:5px 11px;font-size:12px;cursor:pointer}
button:hover{border-color:var(--accent)}
button.on{background:var(--accent-soft);border-color:var(--accent);color:var(--accent-ink)}
button:focus-visible{outline:2px solid var(--accent);outline-offset:2px}
main{padding:16px 22px 30px;display:flex;flex-direction:column;gap:14px}
.card{background:var(--panel);border:1px solid var(--line);border-radius:10px;
      padding:12px 14px 10px;box-shadow:var(--shadow)}
.card h2{margin:0 0 8px;font-size:13px;font-weight:600}
.filters{display:flex;gap:6px;flex-wrap:wrap;align-items:center}
.legend{display:flex;gap:14px;flex-wrap:wrap;font-size:11.5px;color:var(--ink-2);margin:6px 0 8px}
.legend i{display:inline-block;width:11px;height:11px;border-radius:3px;margin-right:5px;
          vertical-align:-1px}
.wrap{position:relative;width:100%;overflow-x:auto}
canvas{display:block;width:100%}
table{border-collapse:collapse;font-size:12.5px;font-variant-numeric:tabular-nums}
td,th{padding:3px 14px 3px 0;text-align:left}
th{color:var(--ink-3);font-weight:500}
.bad{color:var(--fit);font-weight:600}
.ok{color:var(--good)}
.note{font-size:12px;color:var(--ink-2);max-width:78ch}
@media (prefers-reduced-motion:reduce){*{transition:none!important}}
</style>

<header>
  <div>
    <div class="eyebrow">a1975 &middot; AT-TPC &middot; FRIB generic channel 0</div>
    <h1>ion-chamber traces &mdash; how the height is measured</h1>
  </div>
  <div class="spacer"></div>
  <button id="theme">theme</button>
</header>

<main>
  <div class="card">
    <div class="filters" id="filters"></div>
  </div>

  <div class="card">
    <h2 id="ttl">trace</h2>
    <div class="legend">
      <span><i style="background:var(--data)"></i>ADC samples</span>
      <span><i style="background:var(--accent);opacity:.45"></i>amplitude window [1050,1250] &mdash; what icmax uses</span>
      <span><i style="background:var(--ref2)"></i>pulse-count window [800,1500]</span>
      <span><i style="background:var(--fit)"></i>threshold 200 (and half, 100)</span>
    </div>
    <div class="wrap"><canvas id="cv"></canvas></div>
    <div class="filters" style="margin-top:10px">
      <button id="prev">&larr; previous</button>
      <button id="next">next &rarr;</button>
      <span class="mono" id="pos" style="color:var(--ink-3)"></span>
      <span class="spacer"></span>
      <span class="note">arrow keys also work</span>
    </div>
  </div>

  <div class="card">
    <h2>the two reductions, for this trace</h2>
    <table id="tab"></table>
    <p class="note" id="expl" style="margin:10px 0 2px"></p>
  </div>
</main>

<script>
const D = __DATA__;
const T = D.traces;
let filter = 'all', idx = 0;

const CLASSES = ['all', ...Array.from(new Set(T.map(t => t.cls)))];
const view = () => filter === 'all' ? T : T.filter(t => t.cls === filter);

const fb = document.getElementById('filters');
CLASSES.forEach(c => {
  const b = document.createElement('button');
  b.textContent = c === 'all' ? ('all (' + T.length + ')') : (c + ' (' + T.filter(t => t.cls === c).length + ')');
  b.onclick = () => { filter = c; idx = 0; paint(); };
  b.dataset.cls = c;
  fb.appendChild(b);
});

const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();

function draw(tr){
  const cv = document.getElementById('cv');
  const dpr = window.devicePixelRatio || 1;
  const w = cv.clientWidth || 900, h = 340;
  cv.width = w * dpr; cv.height = h * dpr; cv.style.height = h + 'px';
  const x = cv.getContext('2d'); x.setTransform(dpr,0,0,dpr,0,0);
  x.clearRect(0,0,w,h);
  const L = 58, R = 14, TP = 12, B = 34, pw = w-L-R, ph = h-TP-B;
  const n = tr.adc.length;
  const tb0 = D.tbLo, tb1 = D.tbLo + n;
  let amax = Math.max(300, ...tr.adc), amin = Math.min(0, ...tr.adc);
  const pad = (amax-amin)*0.10; amax += pad;
  const px = tb => L + (tb-tb0)/(tb1-tb0)*pw;
  const py = v  => TP + ph - (v-amin)/(amax-amin)*ph;

  // amplitude window (shaded) and pulse window (outlined)
  x.fillStyle = css('--accent-soft');
  x.fillRect(px(D.icTbLo), TP, px(D.icTbHi)-px(D.icTbLo), ph);
  x.strokeStyle = css('--ref2'); x.lineWidth = 1.5; x.setLineDash([5,4]);
  x.strokeRect(px(D.pkTbLo), TP+1, px(D.pkTbHi)-px(D.pkTbLo), ph-2);
  x.setLineDash([]);
  x.strokeStyle = css('--accent'); x.lineWidth = 1;
  [D.icTbLo, D.icTbHi].forEach(t => { x.beginPath(); x.moveTo(px(t),TP); x.lineTo(px(t),TP+ph); x.stroke(); });

  // grid
  x.strokeStyle = css('--grid'); x.lineWidth = 1;
  for (let t = 800; t <= 1500; t += 100){ x.beginPath(); x.moveTo(px(t),TP); x.lineTo(px(t),TP+ph); x.stroke(); }

  // thresholds
  x.strokeStyle = css('--fit'); x.setLineDash([4,3]); x.globalAlpha = .85;
  [D.thr, D.thr/2].forEach(v => { if (v<amax){ x.beginPath(); x.moveTo(L,py(v)); x.lineTo(L+pw,py(v)); x.stroke(); }});
  x.setLineDash([]); x.globalAlpha = 1;

  // the trace
  x.strokeStyle = css('--data'); x.lineWidth = 1.6; x.beginPath();
  for (let i = 0; i < n; ++i){ const X = px(tb0+i), Y = py(tr.adc[i]); i ? x.lineTo(X,Y) : x.moveTo(X,Y); }
  x.stroke();

  // the two maxima
  const mark = (tb, v, col, lab) => {
    x.fillStyle = col; x.beginPath(); x.arc(px(tb), py(v), 4.5, 0, 7); x.fill();
    x.font = '11.5px ui-monospace,SFMono-Regular,"DejaVu Sans Mono",Menlo,monospace';
    x.textAlign = 'left'; x.fillText(lab, Math.min(px(tb)+8, L+pw-120), py(v)-7);
  };
  mark(tr.trueTb, tr.trueMax, css('--data'), 'true max ' + tr.trueMax);
  if (tr.fixMax !== tr.trueMax){
    // where the fixed window's maximum actually sits
    let bt = D.icTbLo, bv = -1e9;
    for (let t = D.icTbLo; t < D.icTbHi; ++t){ const v = tr.adc[t-tb0]; if (v !== undefined && v > bv){ bv = v; bt = t; } }
    mark(bt, bv, css('--fit'), 'icmax ' + tr.fixMax + '  <- what is recorded');
  }

  // axes
  x.strokeStyle = css('--line'); x.lineWidth = 1;
  x.beginPath(); x.moveTo(L,TP); x.lineTo(L,TP+ph); x.lineTo(L+pw,TP+ph); x.stroke();
  x.fillStyle = css('--ink-3');
  x.font = '11px ui-monospace,SFMono-Regular,"DejaVu Sans Mono",Menlo,monospace';
  x.textAlign = 'center';
  for (let t = 800; t <= 1500; t += 100) x.fillText(t, px(t), h-16);
  x.fillText('time bucket', L+pw/2, h-3);
  x.textAlign = 'right';
  x.fillText(Math.round(amax), L-6, TP+10);
  x.fillText('0', L-6, py(0));
  x.save(); x.translate(13, TP+ph/2); x.rotate(-Math.PI/2);
  x.textAlign='center'; x.fillText('ADC', 0, 0); x.restore();
}

function paint(){
  const V = view();
  if (!V.length) return;
  idx = ((idx % V.length) + V.length) % V.length;
  const tr = V[idx];
  document.querySelectorAll('#filters button').forEach(b =>
    b.classList.toggle('on', b.dataset.cls === filter));
  document.getElementById('ttl').textContent =
    D.run + ' \\u00b7 event ' + tr.ev + ' \\u00b7 ' + tr.cls;
  document.getElementById('pos').textContent = (idx+1) + ' / ' + V.length;
  const agree = tr.fixMax === tr.trueMax;
  document.getElementById('tab').innerHTML =
    '<tr><th>true max, over [' + D.pkTbLo + ',' + D.pkTbHi + ']</th><td class="mono">'
      + tr.trueMax + ' ADC at time bucket ' + tr.trueTb + '</td></tr>'
    + '<tr><th>icmax, over [' + D.icTbLo + ',' + D.icTbHi + ']</th><td class="mono '
      + (agree ? 'ok' : 'bad') + '">' + tr.fixMax + ' ADC</td></tr>'
    + '<tr><th>pulses counted, threshold ' + D.thr + '</th><td class="mono">' + tr.np + '</td></tr>'
    + '<tr><th>true max inside the amplitude window?</th><td class="mono '
      + (tr.outside ? 'bad' : 'ok') + '">' + (tr.outside ? 'NO' : 'yes') + '</td></tr>';
  document.getElementById('expl').textContent = tr.outside
    ? 'The pulse arrived outside [' + D.icTbLo + ',' + D.icTbHi + '], so icmax recorded the '
      + 'baseline instead of the pulse height. The pulse counter still saw it, so this event is '
      + 'labelled single-pulse with an amplitude near zero \\u2014 it then fails the beam window '
      + 'and is discarded. This is 8% of events.'
    : (tr.np > 1
       ? 'More than one pulse in [' + D.pkTbLo + ',' + D.pkTbHi + ']. icmax takes the maximum '
         + 'sample in the amplitude window, which is whichever pulse happens to fall there.'
       : 'Single pulse inside the amplitude window: both reductions find the same height, and '
         + 'this event is measured correctly.');
  draw(tr);
}

document.getElementById('prev').onclick = () => { idx--; paint(); };
document.getElementById('next').onclick = () => { idx++; paint(); };
addEventListener('keydown', e => {
  if (e.key === 'ArrowLeft'){ idx--; paint(); }
  if (e.key === 'ArrowRight'){ idx++; paint(); }
});
document.getElementById('theme').onclick = () => {
  const cur = document.documentElement.getAttribute('data-theme')
           || (matchMedia('(prefers-color-scheme:dark)').matches ? 'dark' : 'light');
  document.documentElement.setAttribute('data-theme', cur === 'dark' ? 'light' : 'dark');
  paint();
};
addEventListener('resize', paint);
matchMedia('(prefers-color-scheme:dark)').addEventListener('change', paint);
paint();
</script>
"""


def main():
    if len(sys.argv) < 3:
        raise SystemExit(__doc__)
    data = json.load(open(sys.argv[1], encoding="utf-8"))
    html = PAGE.replace("__TOKENS__", TOKENS).replace(
        "__DATA__", json.dumps(data, separators=(",", ":")))
    open(sys.argv[2], "w", encoding="utf-8").write(html)
    print("wrote %s : %d traces, %.2f MB"
          % (sys.argv[2], len(data["traces"]), len(html.encode()) / 1048576))


if __name__ == "__main__":
    main()
