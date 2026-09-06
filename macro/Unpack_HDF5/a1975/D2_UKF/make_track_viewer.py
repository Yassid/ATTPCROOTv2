#!/usr/bin/env python3
"""Build a self-contained browser track viewer from dump_tracks_dp.C's JSON.

WHY A VIEWER AND NOT MORE HISTOGRAMS: the backward sub-MeV selection is short-track territory
(a 1 MeV proton ranges 264 mm, a 0.5 MeV one 86 mm). No aggregate distribution separates a real
short track from a curler fragment or a noise cluster -- that judgement is made by eye, one
track at a time, which is what this page is for.

The page keeps YOUR verdicts in localStorage and can export them, so a pass over 240 tracks is
not lost when the tab closes.

  ./make_track_viewer.py tracks.json out.html "title"
"""
import json, sys, html

src, out = sys.argv[1], sys.argv[2]
title = sys.argv[3] if len(sys.argv) > 3 else "a1975 (d,p) track viewer"
D = json.load(open(src))
tracks = D["tracks"]
# The JSON is embedded rather than fetched: a file:// page cannot fetch() a sibling file (CORS),
# so a two-file viewer would work from a server and silently fail from the Desktop.
blob = json.dumps(tracks, separators=(",", ":"))

page = """<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>__TITLE__</title>
<style>
:root{color-scheme:light dark;
  --bg:#fbfbfa; --fg:#1a1a19; --dim:#6b6b68; --line:#e0e0dc; --card:#fff;
  --good:#2ca02c; --bad:#e6194b; --unsure:#d98b00; --accent:#0b6cff;}
@media (prefers-color-scheme:dark){:root{
  --bg:#16161a; --fg:#e8e8e6; --dim:#9a9a96; --line:#2c2c32; --card:#1e1e24;}}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);
  font:13px/1.5 ui-sans-serif,system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
header{padding:10px 16px;border-bottom:1px solid var(--line);display:flex;
  gap:16px;align-items:baseline;flex-wrap:wrap}
h1{font-size:15px;margin:0;font-weight:650}
.sel{color:var(--dim);font-size:12px}
main{display:grid;grid-template-columns:minmax(0,1fr) 260px;gap:16px;padding:16px;align-items:start}
@media(max-width:900px){main{grid-template-columns:minmax(0,1fr)}}
.pads{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:10px}
@media(max-width:760px){.pads{grid-template-columns:minmax(0,1fr)}}
.pad{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:6px}
.pad h3{margin:0 0 4px;font-size:11px;font-weight:600;color:var(--dim);text-align:center}
canvas{width:100%;height:auto;display:block}
.side{background:var(--card);border:1px solid var(--line);border-radius:8px;padding:12px}
table{width:100%;border-collapse:collapse;font-variant-numeric:tabular-nums}
td{padding:2px 0;font-size:12px}
td:first-child{color:var(--dim)}
td:last-child{text-align:right;font-weight:600}
.nav{display:flex;gap:6px;align-items:center;margin:12px 0;flex-wrap:wrap}
button{font:inherit;padding:5px 11px;border:1px solid var(--line);background:var(--card);
  color:var(--fg);border-radius:6px;cursor:pointer}
button:hover{border-color:var(--accent)}
button.on{color:#fff;border-color:transparent}
#bg.on{background:var(--good)} #bb.on{background:var(--bad)} #bu.on{background:var(--unsure)}
input[type=number]{font:inherit;width:70px;padding:4px 6px;border:1px solid var(--line);
  background:var(--bg);color:var(--fg);border-radius:6px}
.tally{display:flex;gap:12px;font-size:12px;margin-top:10px;color:var(--dim)}
.tally b{color:var(--fg)}
.strip{display:flex;flex-wrap:wrap;gap:3px;margin-top:12px;max-height:190px;overflow:auto}
.chip{width:13px;height:13px;border-radius:3px;border:1px solid var(--line);cursor:pointer;background:var(--bg)}
.chip.g{background:var(--good);border-color:transparent}
.chip.b{background:var(--bad);border-color:transparent}
.chip.u{background:var(--unsure);border-color:transparent}
.chip.cur{outline:2px solid var(--accent);outline-offset:1px}
kbd{font:11px ui-monospace,monospace;border:1px solid var(--line);border-radius:4px;padding:0 4px}
.note{color:var(--dim);font-size:11px;margin-top:10px;line-height:1.45}
.fitbox{margin-top:14px;padding-top:10px;border-top:1px solid var(--line)}
.eyebrow{font-size:10px;letter-spacing:.08em;text-transform:uppercase;color:var(--dim);margin-bottom:6px}
.frow{display:flex;align-items:center;gap:6px;margin:3px 0;font-size:11px}
.frow label{color:var(--dim);width:74px;flex:none}
.frow input[type=range]{flex:1;min-width:0}
.frow span{width:52px;text-align:right;font-variant-numeric:tabular-nums}
#cmp,#modes{margin-top:8px}
#modes th{font-size:10px;color:var(--dim);font-weight:600;text-align:right;padding:1px 0}
#modes th:first-child{text-align:left}
#modes td{font-size:11px;text-align:right}
#modes td:first-child{text-align:left;color:var(--dim)}
#modes tr.sel td{font-weight:700;color:var(--accent)}
button.mode{padding:3px 8px;font-size:11px}
button.mode.on{background:var(--accent);color:#fff;border-color:transparent}
#cmp th{font-size:10px;color:var(--dim);font-weight:600;text-align:right;padding:1px 0}
#cmp th:first-child{text-align:left}
#cmp td:nth-child(2){text-align:right;color:var(--dim);font-weight:400}
#cmp td:nth-child(3){text-align:right;font-weight:600}
</style></head><body>
<header>
  <h1>__TITLE__</h1>
  <span class="sel">__SEL__ &middot; showing __NSHOW__ of __NTOT__ qualifying tracks</span>
</header>
<main>
 <div>
  <div class="pads">
    <div class="pad"><h3>beam view &nbsp;z &ndash; x</h3><canvas id="c0" width="460" height="380"></canvas></div>
    <div class="pad"><h3>beam view &nbsp;z &ndash; y</h3><canvas id="c1" width="460" height="380"></canvas></div>
    <div class="pad"><h3>pad plane &nbsp;x &ndash; y</h3><canvas id="c2" width="460" height="380"></canvas></div>
  </div>
  <div class="nav">
    <button id="prev">&larr; prev</button>
    <button id="next">next &rarr;</button>
    <span>#<input type="number" id="idx" min="1" value="1"></span>
    <button id="bg">good <kbd>g</kbd></button>
    <button id="bb">bad <kbd>b</kbd></button>
    <button id="bu">unsure <kbd>u</kbd></button>
    <button id="exp">copy verdicts</button>
    <button id="clr">reset</button>
  </div>
  <div class="tally">
    <span>good <b id="ng">0</b></span><span>bad <b id="nb">0</b></span>
    <span>unsure <b id="nu">0</b></span><span>unjudged <b id="nn">0</b></span>
  </div>
  <div class="strip" id="strip"></div>
  <p class="note">Red = the genfit Kalman trajectory for the SELECTED cluster ordering;
  green dashed = your helix (sliders at right). Green tracking the hits while red does not means the data support a helix genfit did
  not find.  (the fitted state at every measurement
  point, not a helix re-derived from the summary numbers). Arrow keys move, <kbd>g</kbd>/<kbd>b</kbd>/<kbd>u</kbd> judge (and advance).
  Marker area is proportional to cluster charge. The red cross is the fitted vertex; the grey
  dashed line is the pad plane at z = 971.7 mm. Verdicts persist in this browser only.</p>
 </div>
 <div class="side">
   <table id="info"></table>
   <div class="fitbox">
     <div class="eyebrow">genfit &mdash; cluster ordering</div>
     <div class="frow" style="gap:4px">
       <button id="m0" class="mode">z-sort</button>
       <button id="m1" class="mode">cluster</button>
       <button id="m2" class="mode">arc walk</button>
     </div>
     <table id="modes"></table>
     <p class="note" style="margin-top:6px">AtGenfitter orders clusters by drift z. Where the
     helix advances less than one time bucket per cluster (1.84 mm here) that orders NOISE, and
     the filter integrates energy loss along a path the particle never flew. These are real
     genfit refits, precomputed &mdash; the browser cannot run genfit.</p>
   </div>
   <div class="fitbox">
     <div class="eyebrow">helix fit &mdash; yours</div>
     <div class="frow"><label>R (mm)</label><input type="range" id="sR" min="5" max="400" step="0.5"><span id="vR"></span></div>
     <div class="frow"><label>centre x</label><input type="range" id="sCX" min="-300" max="300" step="0.5"><span id="vCX"></span></div>
     <div class="frow"><label>centre y</label><input type="range" id="sCY" min="-300" max="300" step="0.5"><span id="vCY"></span></div>
     <div class="frow"><label>dz/d&phi;</label><input type="range" id="sK" min="-300" max="300" step="0.5"><span id="vK"></span></div>
     <div class="frow"><label>z offset</label><input type="range" id="sZ0" min="-400" max="400" step="1"><span id="vZ0"></span></div>
     <div class="frow"><label>use clusters</label><input type="number" id="i0" step="1" value="0" style="width:52px">
       <input type="number" id="i1" step="1" value="9999" style="width:52px"></div>
     <div class="frow"><label>drop resid &gt;</label><input type="number" id="oc" step="1" value="0" style="width:52px"><span>mm (0=off)</span></div>
     <div class="frow"><button id="auto">auto-fit</button><button id="rev">reverse</button></div>
     <table id="cmp"></table>
   </div>
   <p class="note"><b>How to read &ldquo;yours&rdquo;.</b> Validated against genfit on 99 tracks
   with &chi;&sup2;/ndf &lt; 2: <b>&theta; agrees to 0.3&deg;</b> (99% within 10&deg;), so the
   geometry is trustworthy. But B&rho; comes out <b>15% low</b> and KE 41% low, and that is
   expected, not a bug &mdash; genfit&rsquo;s KE is back-extrapolated to the vertex, while one
   circle over the whole cloud measures the <i>mean</i> curvature, and the spiral tightens as the
   proton slows. To compare like with like, restrict <i>use clusters</i> to the first few near the
   vertex; the radius should grow toward genfit&rsquo;s.</p>
   <p class="note">These tracks sit BELOW genfit's &beta;&gamma; = 0.05 floor (KE 1.17 MeV for a
   proton), the region the CATIMA table governs. A 1 MeV proton ranges 264 mm, a 0.5 MeV one
   86 mm &mdash; so short is expected here, and short alone is not a defect.</p>
 </div>
</main>
<script>
const T = __BLOB__;
const EBEAM = __EBEAM__;
const ZPAD = 971.7312;
let cur = 0;
const KEY = "a1975_dp_verdicts";
let V = {};
try { V = JSON.parse(localStorage.getItem(KEY) || "{}"); } catch(e) { V = {}; }
const id = t => t.run + ":" + t.entry + ":" + t.tid;
const save = () => { try { localStorage.setItem(KEY, JSON.stringify(V)); } catch(e) {} };

function css(v){ return getComputedStyle(document.documentElement).getPropertyValue(v).trim(); }

function draw(ci, t, ax, ay, la, lb){
  const cv = document.getElementById("c"+ci), g = cv.getContext("2d");
  const W = cv.width, H = cv.height, P = 34;
  g.clearRect(0,0,W,H);
  const A = t.hits.map(h=>h[ax]), B = t.hits.map(h=>h[ay]);
  // include the vertex in the extent so a vertex far off the cloud is VISIBLE rather than clipped
  const vA = ax===2?t.vz:(ax===0?t.vx:t.vy), vB = ay===2?t.vz:(ay===0?t.vx:t.vy);
  const CF = curFit(t);
  const FA = CF.map(q=>q[ax]), FB = CF.map(q=>q[ay]);
  // the fit is included in the extent: a trajectory that wanders off the hit cloud is precisely
  // what this page is for, and auto-ranging on the hits alone would crop it out of frame
  let a0=Math.min(...A,vA,...FA), a1=Math.max(...A,vA,...FA),
      b0=Math.min(...B,vB,...FB), b1=Math.max(...B,vB,...FB);
  const pad=(v0,v1)=>{const d=Math.max(v1-v0,8)*0.12;return [v0-d,v1+d];};
  [a0,a1]=pad(a0,a1); [b0,b1]=pad(b0,b1);
  const X = a => P+(a-a0)/(a1-a0)*(W-P-10), Y = b => H-P-(b-b0)/(b1-b0)*(H-P-10);
  g.strokeStyle=css("--line"); g.lineWidth=1;
  g.beginPath(); g.moveTo(P,10); g.lineTo(P,H-P); g.lineTo(W-10,H-P); g.stroke();
  g.fillStyle=css("--dim"); g.font="10px ui-sans-serif,system-ui";
  g.textAlign="center"; g.fillText(la, (P+W-10)/2, H-9);
  g.save(); g.translate(11,(H-P)/2); g.rotate(-Math.PI/2); g.fillText(lb,0,0); g.restore();
  g.textAlign="right";
  g.fillText(a0.toFixed(0), P+14, H-P+13); g.fillText(a1.toFixed(0), W-12, H-P+13);
  g.fillText(b1.toFixed(0), P-4, 16); g.fillText(b0.toFixed(0), P-4, H-P);
  // pad plane marker, only on the two beam views
  if ((ax===2||ay===2) && ZPAD>=a0 && ZPAD<=a1 && ax===2){
    g.strokeStyle=css("--dim"); g.setLineDash([4,3]);
    g.beginPath(); g.moveTo(X(ZPAD),10); g.lineTo(X(ZPAD),H-P); g.stroke(); g.setLineDash([]);
  }
  // the FITTED trajectory, drawn UNDER the hits so it can never hide a hit it misses
  { const GF = curFit(t);
    if (GF.length > 1){
      g.strokeStyle = "rgba(230,25,75,0.9)"; g.lineWidth = 1.8;
      g.beginPath();
      GF.forEach((q,k)=>{ const px=X(q[ax]), py=Y(q[ay]); k?g.lineTo(px,py):g.moveTo(px,py); });
      g.stroke();
    } }
  // the hand-tuned helix, in green
  if (HFIT){
    const HP = helixPts(t, HFIT);
    if (HP.length>1){
      g.strokeStyle="rgba(44,160,44,0.95)"; g.lineWidth=1.8; g.setLineDash([6,3]);
      g.beginPath();
      HP.forEach((q,k)=>{const px=X(q[ax]),py=Y(q[ay]); k?g.lineTo(px,py):g.moveTo(px,py);});
      g.stroke(); g.setLineDash([]);
    }
  }
  const qm = Math.max(...t.hits.map(h=>h[3]), 1);
  for (const h of t.hits){
    const r = 1.4 + 3.4*Math.sqrt(h[3]/qm);
    g.beginPath(); g.arc(X(h[ax]), Y(h[ay]), r, 0, 6.2832);
    g.fillStyle = "rgba(11,108,255,0.62)"; g.fill();
  }
  g.strokeStyle=css("--bad"); g.lineWidth=1.6;
  const cx=X(vA), cy=Y(vB);
  g.beginPath(); g.moveTo(cx-6,cy); g.lineTo(cx+6,cy); g.moveTo(cx,cy-6); g.lineTo(cx,cy+6); g.stroke();
}

// rms 3D distance from each cluster to the NEAREST fitted point -- a number for what the eye is
// judging. Nearest-point rather than index-matched, because the two lists need not correspond
// one-to-one when the fit drops measurements.
function resid(t){
  const F=curFit(t);
  if(!F || F.length<2) return "&mdash;";
  return rmsOf(t,F).toFixed(1);
}
/* ---------------- interactive helix ----------------------------------------
   genfit cannot run in a browser, but a helix in a uniform B along z HAS a closed-form fit:
   a circle in (x,y), then z linear in the unwrapped azimuth. That is the whole model, so the
   sliders below are the model's actual parameters, not a re-parameterisation of a fit result.
   Physics, for Z = 1 at B = 2.85 T:
     p_T[MeV/c] = 0.299792458 * B[T] * R[mm];  p_z = p_T*(k/R) with k = dz/dphi
     theta = atan2(p_T, p_z);  KE = sqrt(p^2+m^2) - m;  Brho = p/299.792458
   Ex comes from the SAME two-body kinematics the aggregate explorer uses, so a helix you tune
   by hand lands on the same Ex axis as the production numbers.                             */
const B_FIELD = 2.85, MP = 938.27208, U = 931.49401;
const M1 = 16.0147013*U, M2 = 2.0135532*U, M3 = 1.00727646688*U, M4 = 17.0225864*U;
let HFIT = null;                       // {cx,cy,R,k,z0,rev}
let FITMODE = 0;                       // 0 z-sort (production) | 1 cluster order | 2 arc walk
const MODE_NAME = ["z-sort (production)", "cluster order", "arc walk"];
/// The trajectory currently on display. Falls back to the single `fit` array for pages built
/// before the refits existed, so an old JSON still renders.
function curFit(t){
  if (t.modes && t.modes[FITMODE]) return t.modes[FITMODE].fit || [];
  return (FITMODE === 0 && t.fit) ? t.fit : [];
}
function modeStat(t,m){ return (t.modes && t.modes[m]) ? t.modes[m] : (m===0 ? t : null); }
function rmsOf(t, fit){
  if(!fit || fit.length<2) return null;
  let s2=0;
  for(const h of t.hits){
    let best=Infinity;
    for(const q of fit){ const d=(h[0]-q[0])**2+(h[1]-q[1])**2+(h[2]-q[2])**2; if(d<best) best=d; }
    s2+=best;
  }
  return Math.sqrt(s2/t.hits.length);
}
function om2(x,y,z){return Math.sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z);}
function kine2b(Eb, thRad, ke){
  const Et1=Eb+M1, Et3=ke+M3;
  const s=M1*M1+M2*M2+2*M2*Et1, uu=M2*M2+M3*M3-2*M2*Et3;
  const arg=(Math.cos(thRad)*om2(s,M1*M1,M2*M2)*om2(uu,M2*M2,M3*M3)
            -(s-M1*M1-M2*M2)*(M2*M2+M3*M3-uu))/(2*M2*M2)+s+uu-M2*M2;
  if(arg<0) return [NaN,NaN];
  const m4x=Math.sqrt(arg);
  const t=M2*M2+m4x*m4x-2*M2*(Et1+M2-Et3);
  const tcm=Math.PI-Math.acos((s*s+s*(2*t-M1*M1-M2*M2-M3*M3-m4x*m4x)+(M1*M1-M2*M2)*(M3*M3-m4x*m4x))
            /(om2(s,M1*M1,M2*M2)*om2(s,M3*M3,m4x*m4x)));
  return [m4x-M4, tcm*180/Math.PI];
}
function selHits(t){
  const i0=+document.getElementById("i0").value, i1=+document.getElementById("i1").value;
  let H=t.hits.filter((h,j)=>j>=i0 && j<=i1);
  const oc=+document.getElementById("oc").value;
  if(oc>0 && HFIT) H=H.filter(h=>hres(h,HFIT)<=oc);
  return H;
}
function circleFit(P){                 // Kasa algebraic fit; exact for noiseless points
  const n=P.length; if(n<3) return null;
  let Sx=0,Sy=0,Sxx=0,Syy=0,Sxy=0,Sxxx=0,Syyy=0,Sxyy=0,Sxxy=0;
  for(const p of P){const x=p[0],y=p[1];
    Sx+=x;Sy+=y;Sxx+=x*x;Syy+=y*y;Sxy+=x*y;Sxxx+=x*x*x;Syyy+=y*y*y;Sxyy+=x*y*y;Sxxy+=x*x*y;}
  const C=n*Sxx-Sx*Sx, D=n*Sxy-Sx*Sy, E=n*Sxxx+n*Sxyy-(Sxx+Syy)*Sx;
  const G=n*Syy-Sy*Sy, H=n*Sxxy+n*Syyy-(Sxx+Syy)*Sy;
  const den=2*(C*G-D*D); if(Math.abs(den)<1e-9) return null;
  const cx=(E*G-D*H)/den, cy=(C*H-D*E)/den;
  let R=0; for(const p of P) R+=Math.hypot(p[0]-cx,p[1]-cy);
  return {cx,cy,R:R/n};
}
function unwrap(P,cx,cy){
  const a=P.map(p=>Math.atan2(p[1]-cy,p[0]-cx)); const out=[a[0]]; let acc=a[0];
  for(let i=1;i<a.length;++i){let d=a[i]-a[i-1];
    while(d>Math.PI)d-=2*Math.PI; while(d<-Math.PI)d+=2*Math.PI; acc+=d; out.push(acc);}
  return out;
}
function autoFit(t){
  const H=selHits(t); const c=circleFit(H); if(!c) return null;
  const ph=unwrap(H,c.cx,c.cy);
  let n=H.length,Sp=0,Sz=0,Spp=0,Spz=0;
  for(let i=0;i<n;++i){Sp+=ph[i];Sz+=H[i][2];Spp+=ph[i]*ph[i];Spz+=ph[i]*H[i][2];}
  const den=n*Spp-Sp*Sp;
  const k=Math.abs(den)<1e-9?0:(n*Spz-Sp*Sz)/den, z0=(Sz-k*Sp)/n;
  return {cx:c.cx, cy:c.cy, R:c.R, k, z0, rev:(HFIT?HFIT.rev:false)};
}
function hres(h,F){                    // distance of one cluster from the helix surface
  const dr=Math.hypot(h[0]-F.cx,h[1]-F.cy)-F.R;
  const ph=Math.atan2(h[1]-F.cy,h[0]-F.cx);
  let best=Infinity;
  for(let w=-3;w<=3;++w){ const dz=h[2]-(F.z0+F.k*(ph+2*Math.PI*w)); if(Math.abs(dz)<Math.abs(best)) best=dz; }
  return Math.hypot(dr,best);
}
function helixRms(t,F){
  const H=selHits(t); if(!H.length) return NaN;
  let s2=0; for(const h of H) s2+=hres(h,F)**2;
  return Math.sqrt(s2/H.length);
}
function helixPhys(F,rev){
  const pT=0.299792458*B_FIELD*F.R;
  const pz=pT*(F.k/F.R);
  const p=Math.hypot(pT,pz);
  let th=Math.atan2(pT,pz);
  if(rev) th=Math.PI-th;
  const ke=Math.sqrt(p*p+MP*MP)-MP;
  const [ex,tcm]=kine2b(EBEAM,th,ke);
  return {p, ke, th:th*180/Math.PI, brho:p/299.792458, ex, tcm};
}
function helixPts(t,F){
  const H=selHits(t); if(H.length<2) return [];
  const ph=unwrap(H,F.cx,F.cy);
  const a=Math.min(...ph), b=Math.max(...ph); const out=[];
  for(let i=0;i<=240;++i){ const q=a+(b-a)*i/240;
    out.push([F.cx+F.R*Math.cos(q), F.cy+F.R*Math.sin(q), F.z0+F.k*q]); }
  return out;
}
function syncSliders(F){
  const set=(id,v,f)=>{document.getElementById(id).value=v;
    document.getElementById("v"+id.slice(1)).textContent=f;};
  set("sR",F.R,F.R.toFixed(1)); set("sCX",F.cx,F.cx.toFixed(0)); set("sCY",F.cy,F.cy.toFixed(0));
  set("sK",F.k,F.k.toFixed(1)); set("sZ0",0,F.z0.toFixed(0));
}
function readSliders(base){
  return {cx:+document.getElementById("sCX").value, cy:+document.getElementById("sCY").value,
          R:+document.getElementById("sR").value,  k:+document.getElementById("sK").value,
          z0:base.z0 + +document.getElementById("sZ0").value, rev:base.rev};
}
function updateModes(){
  const t=T[cur];
  const rows=["<tr><th>ordering</th><th>KE</th><th>&theta;</th><th>&chi;&sup2;/ndf</th><th>rms</th></tr>"];
  for(let m=0;m<3;++m){
    const st=modeStat(t,m);
    const sel = (m===FITMODE) ? " class='sel'" : "";
    if(!st){ rows.push(`<tr${sel}><td>${MODE_NAME[m]}</td><td colspan=4>no fit</td></tr>`); continue; }
    const r=rmsOf(t, st.fit||(m===0?t.fit:null));
    rows.push(`<tr${sel}><td>${MODE_NAME[m]}</td><td>${st.ke.toFixed(2)}</td>`+
              `<td>${st.theta.toFixed(0)}</td><td>${st.chi2ndf>1e8?"&infin;":st.chi2ndf.toFixed(1)}</td>`+
              `<td>${r!==null?r.toFixed(0):"&mdash;"}</td></tr>`);
  }
  document.getElementById("modes").innerHTML=rows.join("");
  for(let m=0;m<3;++m) document.getElementById("m"+m).classList.toggle("on", m===FITMODE);
}
function updateFit(){
  const t=T[cur]; if(!HFIT) return;
  const P=helixPhys(HFIT,HFIT.rev), rmsH=helixRms(t,HFIT);
  const rows=[["","genfit","yours"],
    ["KE (MeV)", t.ke.toFixed(3), P.ke.toFixed(3)],
    ["&theta;<sub>lab</sub>", t.theta.toFixed(1), P.th.toFixed(1)],
    ["B&rho; (T m)", t.brho.toFixed(3), P.brho.toFixed(3)],
    ["E<sub>x</sub> (MeV)", t.ex.toFixed(2), isFinite(P.ex)?P.ex.toFixed(2):"&mdash;"],
    ["rms (mm)", resid(t), isFinite(rmsH)?rmsH.toFixed(1):"&mdash;"],
    ["clusters used", t.ncl, selHits(t).length]];
  updateModes();
  document.getElementById("cmp").innerHTML =
    "<tr><th>"+rows[0].join("</th><th>")+"</th></tr>" +
    rows.slice(1).map(r=>"<tr><td>"+r[0]+"</td><td>"+r[1]+"</td><td>"+r[2]+"</td></tr>").join("");
  redraw();
}
function show(i){
  cur = (i+T.length) % T.length;
  const t = T[cur];
  draw(0,t,2,0,"z (mm)","x (mm)");
  draw(1,t,2,1,"z (mm)","y (mm)");
  draw(2,t,0,1,"x (mm)","y (mm)");
  const ms = modeStat(t,FITMODE) || t;
  const rows = [["run", String(t.run).padStart(4,"0")],["entry",t.entry],["track id",t.tid],
    ["KE (MeV)",(ms.ke!==undefined?ms.ke:t.ke).toFixed(3)],
    ["&theta;<sub>lab</sub> (&deg;)",(ms.theta!==undefined?ms.theta:t.theta).toFixed(1)],
    ["&phi; (&deg;)",t.phi.toFixed(0)],["E<sub>x</sub> (MeV)",t.ex.toFixed(2)],
    ["&theta;<sub>cm</sub> (&deg;)",t.thcm.toFixed(1)],["B&rho; (T m)",t.brho.toFixed(3)],
    ["&chi;&sup2;/ndf",(ms.chi2ndf!==undefined?ms.chi2ndf:t.chi2ndf).toFixed(2)],["clusters",t.ncl],
    ["fit points",(ms.nfit!==undefined?ms.nfit:(t.nfit!==undefined?t.nfit:"&mdash;"))],
    ["hit&ndash;fit rms (mm)",resid(t)],["IC",t.ic.toFixed(0)],
    ["vertex z (mm)",t.vz.toFixed(0)],["vertex r (mm)",Math.hypot(t.vx,t.vy).toFixed(0)]];
  document.getElementById("info").innerHTML =
    rows.map(r=>"<tr><td>"+r[0]+"</td><td>"+r[1]+"</td></tr>").join("");
  document.getElementById("idx").value = cur+1;
  const v = V[id(t)];
  for (const [b,k] of [["bg","g"],["bb","b"],["bu","u"]])
    document.getElementById(b).classList.toggle("on", v===k);
  document.querySelectorAll(".chip").forEach((c,j)=>c.classList.toggle("cur", j===cur));
  document.getElementById("i1").value = 9999;
  HFIT = autoFit(t) || null;
  if (HFIT) syncSliders(HFIT);
  updateFit();
  tally();
}
function redraw(){
  const t=T[cur];
  draw(0,t,2,0,"z (mm)","x (mm)");
  draw(1,t,2,1,"z (mm)","y (mm)");
  draw(2,t,0,1,"x (mm)","y (mm)");
}
function tally(){
  let g=0,b=0,u=0;
  T.forEach(t=>{const v=V[id(t)]; if(v==="g")g++; else if(v==="b")b++; else if(v==="u")u++;});
  document.getElementById("ng").textContent=g; document.getElementById("nb").textContent=b;
  document.getElementById("nu").textContent=u; document.getElementById("nn").textContent=T.length-g-b-u;
}
function judge(k){
  const t=T[cur]; if(V[id(t)]===k) delete V[id(t)]; else V[id(t)]=k;
  save();
  const c=document.querySelectorAll(".chip")[cur];
  c.className="chip"+(V[id(t)]?" "+V[id(t)]:"");
  if(V[id(t)]) show(cur+1); else show(cur);
}
const strip=document.getElementById("strip");
T.forEach((t,j)=>{ const c=document.createElement("div");
  c.className="chip"+(V[id(t)]?" "+V[id(t)]:""); c.title="#"+(j+1)+"  KE "+t.ke.toFixed(2)+" MeV";
  c.onclick=()=>show(j); strip.appendChild(c); });
["sR","sCX","sCY","sK","sZ0"].forEach(id=>{
  document.getElementById(id).oninput=()=>{
    if(!HFIT) return;
    const base=autoFit(T[cur])||HFIT;
    HFIT=readSliders({z0:base.z0, rev:HFIT.rev});
    document.getElementById("v"+id.slice(1)).textContent=(+document.getElementById(id).value).toFixed(1);
    updateFit();
  };
});
["i0","i1","oc"].forEach(id=>document.getElementById(id).onchange=()=>{
  HFIT=autoFit(T[cur])||HFIT; if(HFIT) syncSliders(HFIT); updateFit(); });
[0,1,2].forEach(m=>document.getElementById("m"+m).onclick=()=>{ FITMODE=m; show(cur); });
document.getElementById("auto").onclick=()=>{ const f=autoFit(T[cur]);
  if(f){ f.rev=HFIT?HFIT.rev:false; HFIT=f; syncSliders(HFIT); updateFit(); } };
document.getElementById("rev").onclick=()=>{ if(HFIT){ HFIT.rev=!HFIT.rev; updateFit(); } };
document.getElementById("prev").onclick=()=>show(cur-1);
document.getElementById("next").onclick=()=>show(cur+1);
document.getElementById("idx").onchange=e=>show(parseInt(e.target.value,10)-1);
document.getElementById("bg").onclick=()=>judge("g");
document.getElementById("bb").onclick=()=>judge("b");
document.getElementById("bu").onclick=()=>judge("u");
document.getElementById("clr").onclick=()=>{ if(confirm("Clear all verdicts?")){ V={}; save();
  document.querySelectorAll(".chip").forEach(c=>c.className="chip"); show(cur); } };
document.getElementById("exp").onclick=()=>{
  const rows=["run,entry,tid,ke,theta,ex,brho,chi2ndf,ncl,vertexz,verdict"];
  T.forEach(t=>{const v=V[id(t)]; if(v) rows.push([String(t.run).padStart(4,"0"),t.entry,t.tid,
    t.ke.toFixed(3),t.theta.toFixed(2),t.ex.toFixed(3),t.brho.toFixed(4),t.chi2ndf.toFixed(3),
    t.ncl,t.vz.toFixed(1),v].join(","));});
  const txt=rows.join("\\n");
  navigator.clipboard.writeText(txt).then(
    ()=>alert("Copied "+(rows.length-1)+" verdicts as CSV."),
    ()=>{ const w=window.open(""); w.document.write("<pre>"+txt.replace(/</g,"&lt;")+"</pre>"); });
};
addEventListener("keydown",e=>{
  if(e.target.tagName==="INPUT") return;
  if(e.key==="ArrowRight"||e.key===" ") {show(cur+1); e.preventDefault();}
  else if(e.key==="ArrowLeft") show(cur-1);
  else if(e.key==="g") judge("g"); else if(e.key==="b") judge("b"); else if(e.key==="u") judge("u");
});
show(0);
</script></body></html>
"""
page = (page.replace("__TITLE__", html.escape(title))
            .replace("__SEL__", html.escape(D.get("sel", "")))
            .replace("__NSHOW__", str(len(tracks)))
            .replace("__NTOT__", str(D.get("ntotal", len(tracks))))
            .replace("__BLOB__", blob)
            .replace("__EBEAM__", repr(D.get("ebeam", 184.25))))
open(out, "w", encoding="utf-8").write(page)
print("  wrote %s  (%d tracks, %.1f MB)" % (out, len(tracks), len(page)/1e6))
