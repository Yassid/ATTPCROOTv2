/// @file make_gate_editor.C
/// @brief Bake the ungated Spyral PID landscape into a self-contained browser gate editor.
///
/// The X11 editor (pp/pid_gate_edit.C) reads the FIT output, which only holds tracks that already
/// passed the gate being refined -- you cannot see what you are cutting away. This one is built on
/// the ungated landscape from pid/pid_plane_cache.C, so the proton / deuteron / triton bands are
/// all visible while the polygon is drawn.
///
/// The page bakes a 2D HISTOGRAM, not the raw tracks: ~1M tracks would be a ~30 MB JSON, while a
/// sparse histogram is a few hundred kB and counts inside a polygon just as well (the count is a
/// sum over bins whose CENTRE is inside, so it is accurate to the bin size, which is far finer
/// than any sensible gate edge).
///
/// Three nCluster stacks are baked (>=0, >=15, >=25) so track quality can be tightened in-page:
/// the low-quality tracks are what smear the bands together.
///
/// Output JSON is exactly the AtCut2D schema (name / xaxis / yaxis / vertices), so it can be fed
/// straight to AtGenfitter::SetPIDGate or AtParticleID::LoadJSON.
///
///   root -b -q 'pid/make_gate_editor.C("/tmp/pidplane/pid_plane.root","/home/yassid/pid_gate.html")'

#include <cstdio>
#include <string>
#include <vector>

static std::string readAll(const char *p)
{
   FILE *f = fopen(p, "rb");
   if (!f) return "";
   std::string s;
   char buf[65536]; size_t n;
   while ((n = fread(buf, 1, sizeof buf, f)) > 0) s.append(buf, n);
   fclose(f);
   return s;
}

void make_gate_editor(TString cache = "", TString outHtml = "/home/yassid/a1975_pid_gate.html",
                      TString tag = "16C+p  (a1975)", double xMax = 40.0, double yMax = 2.5,
                      int nx = 400, int ny = 250, TString gatesCSV = "pid/deuteron_band.json,pid/proton_band.json")
{
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   if (cache.IsNull()) { printf("give the pid_plane.root cache\n"); return; }
   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", cache.Data()); return; }
   TTree *t = (TTree *)f->Get("pl");
   if (!t) { printf("no tree `pl` in %s\n", cache.Data()); return; }

   float sqrtdedx, brho, polar; int nclusters;
   t->SetBranchAddress("sqrtdedx",&sqrtdedx); t->SetBranchAddress("brho",&brho);
   t->SetBranchAddress("polar",&polar);       t->SetBranchAddress("nclusters",&nclusters);

   const int NC = 3;
   const int cut[NC] = {0, 15, 25};
   std::vector<std::vector<int>> H(NC, std::vector<int>((size_t)nx*ny, 0));
   std::vector<long> tot(NC, 0);
   const Long64_t N = t->GetEntries();
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      if (sqrtdedx < 0 || sqrtdedx >= xMax || brho < 0 || brho >= yMax) continue;
      int ix = (int)(sqrtdedx/xMax*nx), iy = (int)(brho/yMax*ny);
      if (ix < 0 || ix >= nx || iy < 0 || iy >= ny) continue;
      for (int c = 0; c < NC; ++c)
         if (nclusters >= cut[c]) { H[c][(size_t)iy*nx + ix]++; tot[c]++; }
   }
   f->Close();
   printf("landscape: %lld tracks binned; totals", N);
   for (int c = 0; c < NC; ++c) printf("  nclus>=%d: %ld", cut[c], tot[c]);
   printf("\n");

   // ---- existing gates, passed through verbatim for overlay ----
   TString gatesJson = "{";
   {
      TObjArray *toks = gatesCSV.Tokenize(",");
      bool first = true;
      for (int i = 0; i < toks->GetEntries(); ++i) {
         TString p = ((TObjString *)toks->At(i))->GetString().Strip(TString::kBoth);
         TString full = p.BeginsWith("/") ? p : here + "/../" + p;
         std::string body = readAll(full.Data());
         if (body.empty()) { printf("  (gate not found, skipped: %s)\n", full.Data()); continue; }
         gatesJson += TString::Format("%s\"%s\":%s", first ? "" : ",",
                                      gSystem->BaseName(p.Data()), body.c_str());
         first = false;
         printf("  overlay gate: %s\n", full.Data());
      }
   }
   gatesJson += "}";

   // ---- sparse histogram payload ----
   TString data = "[";
   for (int c = 0; c < NC; ++c) {
      data += (c ? ",[" : "[");
      bool first = true;
      for (int iy = 0; iy < ny; ++iy)
         for (int ix = 0; ix < nx; ++ix) {
            int v = H[c][(size_t)iy*nx + ix];
            if (!v) continue;
            data += TString::Format("%s%d,%d,%d", first ? "" : ",", ix, iy, v);
            first = false;
         }
      data += "]";
   }
   data += "]";

   FILE *o = fopen(outHtml.Data(), "w");
   if (!o) { printf("cannot write %s\n", outHtml.Data()); return; }
   fprintf(o, "%s", R"HTML(<!doctype html><html><head><meta charset="utf-8">
<title>AT-TPC PID gate editor</title><style>
:root{--bg:#12141a;--fg:#e8eaf0;--dim:#8b93a7;--acc:#e8933a;--pan:#1b1e26;--line:#2a2f3a}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);
font:13px/1.5 ui-monospace,SFMono-Regular,Menlo,monospace}
header{padding:14px 20px;border-bottom:1px solid var(--line)}
h1{margin:0;font-size:17px;letter-spacing:.3px}
.sub{color:var(--dim);font-size:12px;margin-top:3px}
.wrap{display:flex;gap:16px;padding:16px;align-items:flex-start;flex-wrap:wrap}
.panel{background:var(--pan);border:1px solid var(--line);border-radius:8px;padding:14px;min-width:250px}
.panel h2{font-size:11px;text-transform:uppercase;letter-spacing:.12em;color:var(--dim);
margin:0 0 10px;font-weight:600}
label{display:block;margin:9px 0 3px;color:var(--dim);font-size:12px}
input[type=number],select,textarea{width:100%;background:#0d0f14;color:var(--fg);
border:1px solid var(--line);border-radius:5px;padding:6px 8px;font:inherit}
button{background:var(--acc);color:#12141a;border:0;border-radius:5px;padding:7px 11px;
font:inherit;font-weight:700;cursor:pointer;margin:3px 3px 0 0}
button.sec{background:#2a2f3a;color:var(--fg);font-weight:500}
canvas{border:1px solid var(--line);border-radius:8px;cursor:crosshair;background:#0d0f14}
.stat{display:flex;gap:14px;flex-wrap:wrap;margin-top:10px}
.stat div{background:#0d0f14;border:1px solid var(--line);border-radius:6px;padding:8px 12px}
.stat b{display:block;font-size:17px;color:var(--acc)}
.stat span{font-size:11px;color:var(--dim);text-transform:uppercase;letter-spacing:.08em}
.hint{color:var(--dim);font-size:11.5px;margin-top:9px;line-height:1.6}
textarea{height:150px;font-size:11px}
</style></head><body>
<header><h1>AT-TPC PID gate editor</h1><div class="sub" id="sub"></div></header>
<div class="wrap">
 <div class="panel" style="flex:0 0 260px">
  <h2>display</h2>
  <label>track quality</label>
  <select id="q"><option value="0">all tracks</option>
  <option value="1" selected>nClusters &ge; 15</option><option value="2">nClusters &ge; 25</option></select>
  <label><input type="checkbox" id="logz" checked> log colour scale</label>
  <label>saturate at (counts)</label><input type="number" id="sat" value="0" step="1">
  <div class="hint">0 = auto (99.5th percentile)</div>
  <label>overlay existing gate</label><select id="ov"></select>
  <h2 style="margin-top:16px">gate</h2>
  <button id="undo" class="sec">undo point</button>
  <button id="clear" class="sec">clear</button>
  <button id="close">close polygon</button>
  <button id="seed" class="sec">start from overlay</button>
  <label>name</label><input type="text" id="gname" value="deuteron_tight"
   style="width:100%;background:#0d0f14;color:var(--fg);border:1px solid var(--line);border-radius:5px;padding:6px 8px;font:inherit">
  <label>Z, A (optional)</label>
  <div style="display:flex;gap:6px"><input type="number" id="gz" value="1"><input type="number" id="ga" value="2"></div>
  <button id="exp">export JSON</button>
  <div class="hint">Click the plot to add vertices. The count updates live.
  Export writes the AtCut2D schema, ready for AtGenfitter::SetPIDGate.</div>
 </div>
 <div>
  <canvas id="cv" width="900" height="620"></canvas>
  <div class="stat">
   <div><span>inside gate</span><b id="nin">0</b></div>
   <div><span>total shown</span><b id="ntot">0</b></div>
   <div><span>fraction</span><b id="frac">0 %</b></div>
   <div><span>vertices</span><b id="nv">0</b></div>
  </div>
 </div>
 <div class="panel" style="flex:1 1 320px">
  <h2>json</h2><textarea id="out" readonly></textarea>
  <button id="copy" class="sec">copy to clipboard</button>
  <button id="dl" class="sec">download</button>
  <div class="hint">Save it into <code>macro/Unpack_HDF5/a1975/UKF/pid/</code> and pass the path
  as the <code>pidGate</code> argument of the fit macro.</div>
 </div>
</div>
<script>
const XMAX=__XMAX__, YMAX=__YMAX__, NX=__NX__, NY=__NY__, TAG="__TAG__";
const TOT=__TOT__, SPARSE=__DATA__, GATES=__GATES__;
document.getElementById('sub').textContent =
  TAG+" · Spyral PID: √dE/dx vs Bρ · ungated pattern tracks";
// rebuild dense grids from the sparse [ix,iy,v,...] triplets
const GRID=SPARSE.map(a=>{const g=new Int32Array(NX*NY);
  for(let i=0;i<a.length;i+=3) g[a[i+1]*NX+a[i]]=a[i+2]; return g;});
const cv=document.getElementById('cv'), cx=cv.getContext('2d');
const L=62, R=14, T=12, B=44;                 // plot margins
const PW=cv.width-L-R, PH=cv.height-T-B;
let poly=[], closed=false;
const ovSel=document.getElementById('ov');
ovSel.innerHTML='<option value="">none</option>'+
  Object.keys(GATES).map(k=>`<option value="${k}">${k}</option>`).join('');

const gx=v=>L+v/XMAX*PW, gy=v=>T+PH-v/YMAX*PH;         // data -> pixel
const ix2v=i=>(i+0.5)*XMAX/NX, iy2v=i=>(i+0.5)*YMAX/NY;
const px2x=p=>(p-L)/PW*XMAX, py2y=p=>(T+PH-p)/PH*YMAX;

function inPoly(x,y,P){let c=false;
  for(let i=0,j=P.length-1;i<P.length;j=i++){
    const a=P[i],b=P[j];
    if(((a[1]>y)!=(b[1]>y)) && (x<(b[0]-a[0])*(y-a[1])/(b[1]-a[1])+a[0])) c=!c;}
  return c;}

function draw(){
  const q=+document.getElementById('q').value, g=GRID[q];
  const logz=document.getElementById('logz').checked;
  let sat=+document.getElementById('sat').value;
  if(!(sat>0)){ const nz=[]; for(let i=0;i<g.length;i++) if(g[i]) nz.push(g[i]);
    nz.sort((a,b)=>a-b); sat=nz.length?nz[Math.floor(nz.length*0.995)]:1; }
  cx.clearRect(0,0,cv.width,cv.height);
  cx.fillStyle='#0d0f14'; cx.fillRect(L,T,PW,PH);
  const bw=PW/NX, bh=PH/NY;
  for(let iy=0;iy<NY;iy++)for(let ix=0;ix<NX;ix++){
    const v=g[iy*NX+ix]; if(!v) continue;
    let f=logz?Math.log10(1+v)/Math.log10(1+sat):v/sat; f=Math.max(0,Math.min(1,f));
    // viridis-ish ramp
    const r=Math.round(255*Math.min(1,Math.max(0,-0.4+1.8*f-0.3*f*f)));
    const gg=Math.round(255*Math.min(1,Math.max(0,0.05+1.15*f-0.35*f*f)));
    const b2=Math.round(255*Math.min(1,Math.max(0,0.35+0.9*f-1.4*f*f)));
    cx.fillStyle=`rgb(${r},${gg},${b2})`;
    cx.fillRect(L+ix*bw, T+PH-(iy+1)*bh, Math.ceil(bw), Math.ceil(bh));
  }
  // axes
  cx.strokeStyle='#2a2f3a'; cx.fillStyle='#8b93a7'; cx.lineWidth=1;
  cx.strokeRect(L,T,PW,PH); cx.font='11px ui-monospace,monospace';
  cx.textAlign='center';
  for(let v=0;v<=XMAX;v+=XMAX/8){cx.beginPath();cx.moveTo(gx(v),T+PH);cx.lineTo(gx(v),T+PH+4);
    cx.stroke();cx.fillText(v.toFixed(0),gx(v),T+PH+17);}
  cx.fillText('√dE/dx  [arb]', L+PW/2, T+PH+36);
  cx.textAlign='right';
  for(let v=0;v<=YMAX;v+=YMAX/5){cx.beginPath();cx.moveTo(L,gy(v));cx.lineTo(L-4,gy(v));
    cx.stroke();cx.fillText(v.toFixed(2),L-8,gy(v)+4);}
  cx.save();cx.translate(14,T+PH/2);cx.rotate(-Math.PI/2);cx.textAlign='center';
  cx.fillText('Bρ  [T·m]',0,0);cx.restore();
  // overlay gate
  const ov=ovSel.value;
  if(ov&&GATES[ov]){const V=GATES[ov].vertices;
    cx.strokeStyle='#6ec7ff';cx.lineWidth=1.5;cx.setLineDash([5,3]);cx.beginPath();
    V.forEach((p,i)=>i?cx.lineTo(gx(p[0]),gy(p[1])):cx.moveTo(gx(p[0]),gy(p[1])));
    cx.closePath();cx.stroke();cx.setLineDash([]);}
  // the gate being drawn
  if(poly.length){cx.strokeStyle='#7CFF9B';cx.lineWidth=2;cx.beginPath();
    poly.forEach((p,i)=>i?cx.lineTo(gx(p[0]),gy(p[1])):cx.moveTo(gx(p[0]),gy(p[1])));
    if(closed)cx.closePath(); cx.stroke();
    cx.fillStyle='#7CFF9B'; poly.forEach(p=>cx.fillRect(gx(p[0])-2.5,gy(p[1])-2.5,5,5));}
  count(g);
}

function count(g){
  let nin=0,ntot=0;
  for(let iy=0;iy<NY;iy++)for(let ix=0;ix<NX;ix++){
    const v=g[iy*NX+ix]; if(!v) continue; ntot+=v;
    if(poly.length>2 && inPoly(ix2v(ix),iy2v(iy),poly)) nin+=v;}
  document.getElementById('nin').textContent=nin.toLocaleString();
  document.getElementById('ntot').textContent=ntot.toLocaleString();
  document.getElementById('frac').textContent=ntot?(100*nin/ntot).toFixed(2)+' %':'0 %';
  document.getElementById('nv').textContent=poly.length;
}

cv.addEventListener('click',e=>{
  const r=cv.getBoundingClientRect();
  const x=px2x((e.clientX-r.left)*cv.width/r.width), y=py2y((e.clientY-r.top)*cv.height/r.height);
  if(x<0||x>XMAX||y<0||y>YMAX) return;
  if(closed){poly=[];closed=false;}
  poly.push([+x.toFixed(3),+y.toFixed(4)]); draw();});
document.getElementById('undo').onclick=()=>{poly.pop();closed=false;draw();};
document.getElementById('clear').onclick=()=>{poly=[];closed=false;draw();};
document.getElementById('close').onclick=()=>{if(poly.length>2)closed=true;draw();};
document.getElementById('seed').onclick=()=>{const ov=ovSel.value;
  if(ov&&GATES[ov]){poly=GATES[ov].vertices.map(p=>[p[0],p[1]]);closed=true;draw();}};
['q','logz','sat','ov'].forEach(id=>document.getElementById(id).addEventListener('input',draw));

function json(){
  const o={name:document.getElementById('gname').value,
           xaxis:"sqrtdedx", yaxis:"brho", vertices:poly};
  const z=+document.getElementById('gz').value, a=+document.getElementById('ga').value;
  if(z) o.Z=z; if(a) o.A=a;
  return JSON.stringify(o,null,1);}
document.getElementById('exp').onclick=()=>{document.getElementById('out').value=json();};
document.getElementById('copy').onclick=()=>{const t=document.getElementById('out');
  t.value=t.value||json(); t.select(); document.execCommand('copy');};
document.getElementById('dl').onclick=()=>{const b=new Blob([json()],{type:'application/json'});
  const u=URL.createObjectURL(b), a=document.createElement('a');
  a.href=u; a.download=document.getElementById('gname').value+'.json'; a.click();
  URL.revokeObjectURL(u);};
draw();
</script></body></html>)HTML");
   fclose(o);

   // splice the payload in (kept out of the raw string so the HTML stays readable)
   std::string html = readAll(outHtml.Data());
   auto rep = [&](const std::string &k, const std::string &v) {
      size_t p;
      while ((p = html.find(k)) != std::string::npos) html.replace(p, k.size(), v);
   };
   rep("__XMAX__", TString::Format("%g", xMax).Data());
   rep("__YMAX__", TString::Format("%g", yMax).Data());
   rep("__NX__", TString::Format("%d", nx).Data());
   rep("__NY__", TString::Format("%d", ny).Data());
   rep("__TAG__", tag.Data());
   rep("__TOT__", TString::Format("%ld", tot[1]).Data());
   rep("__DATA__", data.Data());
   rep("__GATES__", gatesJson.Data());
   o = fopen(outHtml.Data(), "w");
   fwrite(html.data(), 1, html.size(), o);
   fclose(o);
   printf("wrote %s  (%.1f MB)\n", outHtml.Data(), html.size()/1e6);
}
