/// @file make_gate_editor.C
/// @brief Bake the ungated Spyral PID landscape into a self-contained browser gate editor, with a
///        LIVE Brho-vs-angle panel for whatever the polygon currently selects.
///
/// The X11 editor (pp/pid_gate_edit.C) reads the FIT output, which only holds tracks that already
/// passed the gate being refined -- you cannot see what you are cutting away. This one is built on
/// the ungated landscape from pid/pid_plane_cache.C, so the proton / deuteron / triton bands are
/// all visible while the polygon is drawn.
///
/// The second panel is the point: a gate in (sqrt(dE/dx), Brho) is only trustworthy if the tracks
/// it selects land on the right KINEMATIC LOCUS, and that is a statement about Brho vs polar angle.
/// Drawing the two side by side turns gate-making from guesswork into a check -- select the wrong
/// band and the selected cloud visibly leaves its curve.
///
/// Data layout: one quantised 4-byte record per track (sqrtdedx, brho, polar, nClusters), base64'd.
/// A 2D histogram would have been smaller but cannot answer "what is the ANGLE distribution of the
/// tracks inside this polygon", which needs the per-track correlation. 520k tracks cost ~2.8 MB.
///
/// Output JSON is exactly the AtCut2D schema (name / xaxis / yaxis / vertices), so it can be fed
/// straight to AtGenfitter::SetPIDGate or AtParticleID::LoadJSON.
///
///   root -b -q 'pid/make_gate_editor.C("/tmp/pidplane_d2/pid_plane_d2.root","/home/yassid/g.html")'

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

static std::string b64(const std::vector<unsigned char> &v)
{
   static const char *T = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
   std::string o;
   o.reserve((v.size() + 2)/3*4);
   for (size_t i = 0; i < v.size(); i += 3) {
      unsigned x = v[i] << 16;
      if (i+1 < v.size()) x |= v[i+1] << 8;
      if (i+2 < v.size()) x |= v[i+2];
      o += T[(x>>18)&63]; o += T[(x>>12)&63];
      o += (i+1 < v.size()) ? T[(x>>6)&63] : '=';
      o += (i+2 < v.size()) ? T[x&63] : '=';
   }
   return o;
}

static double om_(double x, double y, double z) { return sqrt(x*x+y*y+z*z-2*x*y-2*y*z-2*x*z); }
static double exOf_(double m1, double m2, double m3, double m4, double K, double th, double Ke)
{
   double E1 = K + m1, E3 = Ke + m3;
   double s = m1*m1+m2*m2+2*m2*E1, u = m2*m2+m3*m3-2*m2*E3;
   double a = (cos(th)*om_(s,m1*m1,m2*m2)*om_(u,m2*m2,m3*m3)-(s-m1*m1-m2*m2)*(m2*m2+m3*m3-u))/(2*m2*m2)
              + s + u - m2*m2;
   return a < 0 ? std::nan("") : sqrt(a) - m4;
}
/// All Brho values with Ex(KE)=exT at this theta_lab.
///
/// Ex(KE) at fixed angle is NOT monotonic once the ejectile is heavy: for 16C(d,t)15C it rises then
/// falls, so a given Ex has TWO solutions at forward angles and NONE past the channel's maximum lab
/// angle (~40 deg here). A single bisection -- correct for a light ejectile like (p,p') -- silently
/// dropped the entire triton locus. So scan for sign changes and bisect each bracket separately;
/// an empty result is a real statement that the channel cannot reach that angle.
static std::vector<double> brhoAt_(double m1, double m2, double m3, double m4, double E, double th,
                                   double exT, int Z)
{
   std::vector<double> out;
   const int NS = 600;
   const double keMax = 300.0;
   double prevKe = 0, prevF = std::nan("");
   for (int i = 0; i <= NS; ++i) {
      double ke = 0.05 + (keMax - 0.05)*i/NS;
      double f = exOf_(m1,m2,m3,m4,E,th,ke);
      if (!std::isnan(f) && !std::isnan(prevF) && (prevF-exT)*(f-exT) <= 0) {
         double lo = prevKe, hi = ke, flo = prevF;
         for (int k = 0; k < 60; ++k) {
            double mid = 0.5*(lo+hi), fm = exOf_(m1,m2,m3,m4,E,th,mid);
            if (std::isnan(fm)) { hi = mid; continue; }
            if ((flo-exT)*(fm-exT) <= 0) hi = mid; else { lo = mid; flo = fm; }
         }
         double k2 = 0.5*(lo+hi), p = sqrt(k2*k2 + 2*k2*m3);
         out.push_back(p/(299.792458*Z));
      }
      prevKe = ke; prevF = f;
   }
   return out;
}

/// lociSpec: "label|mEjectAmu|mResidAmu|Z;..." -- one kinematic curve per entry, all at Ex=0.
/// Fields are '|'-separated, NOT comma: the labels themselves contain commas, e.g. "(d,t)15C".
void make_gate_editor(TString cache = "", TString outHtml = "/home/yassid/a1975_pid_gate.html",
                      TString tag = "16C+p  (a1975)", double xMax = 40.0, double yMax = 2.5,
                      // 256, not more: the payload is quantised to 256 levels per axis, so a finer
                      // binning leaves periodic empty columns that render as vertical stripes
                      int nx = 256, int ny = 256,
                      TString gatesCSV = "pid/deuteron_band.json,pid/proton_band.json",
                      double Ebeam = 192.0, double mBeamAmu = 16.0147013, double mTargAmu = 2.0135532,
                      TString lociSpec = "(d,p)17C|1.00782503|17.0225864|1;"
                                         "(d,d)16C|2.0135532|16.0147013|1;"
                                         "(d,t)15C|3.01550072|15.0105993|1")
{
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   if (cache.IsNull()) { printf("give the pid_plane*.root cache\n"); return; }
   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", cache.Data()); return; }
   TTree *t = (TTree *)f->Get("pl");
   if (!t) { printf("no tree `pl` in %s\n", cache.Data()); return; }

   float sqrtdedx, brho, polar; int nclusters;
   t->SetBranchAddress("sqrtdedx",&sqrtdedx); t->SetBranchAddress("brho",&brho);
   t->SetBranchAddress("polar",&polar);       t->SetBranchAddress("nclusters",&nclusters);

   std::vector<unsigned char> rec;
   const Long64_t N = t->GetEntries();
   rec.reserve(4*N);
   long kept = 0;
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      if (sqrtdedx < 0 || sqrtdedx >= xMax || brho < 0 || brho >= yMax) continue;
      double pol = polar; if (pol < 0) pol = 0; if (pol > 180) pol = 180;
      rec.push_back((unsigned char)(sqrtdedx/xMax*255.0));
      rec.push_back((unsigned char)(brho/yMax*255.0));
      rec.push_back((unsigned char)(pol/180.0*255.0));
      rec.push_back((unsigned char)(nclusters > 255 ? 255 : (nclusters < 0 ? 0 : nclusters)));
      ++kept;
   }
   f->Close();
   printf("landscape: %lld tracks -> %ld inside the axes (%.1f MB payload)\n",
          N, kept, 4.0*kept/3*4/3e6);

   // kinematic loci, computed here so the page needs no nuclear masses
   const double u = 931.49401;
   TString loci = "[";
   {
      TObjArray *ch = lociSpec.Tokenize(";");
      for (int i = 0; i < ch->GetEntries(); ++i) {
         TObjArray *p = ((TObjString *)ch->At(i))->GetString().Tokenize("|");
         if (p->GetEntries() < 4) continue;
         TString lab = ((TObjString *)p->At(0))->GetString();
         double m3 = ((TObjString *)p->At(1))->GetString().Atof()*u;
         double m4 = ((TObjString *)p->At(2))->GetString().Atof()*u;
         int Z = ((TObjString *)p->At(3))->GetString().Atoi();
         // up to two branches (low-KE and high-KE roots); emit each as its own polyline so the
         // page does not draw a spurious segment joining them
         TString br[2] = {"", ""};
         double aMax = 0;
         for (double a = 1; a <= 179; a += 1) {
            auto sols = brhoAt_(mBeamAmu*u, mTargAmu*u, m3, m4, Ebeam, a*TMath::DegToRad(), 0.0, Z);
            if (!sols.empty()) aMax = a;
            std::sort(sols.begin(), sols.end());
            for (size_t k = 0; k < sols.size() && k < 2; ++k) {
               double b = sols[k];
               if (b <= 0 || b > yMax*1.5) continue;
               br[k] += TString::Format("%s[%.1f,%.4f]", br[k].Length() ? "," : "", a, b);
            }
         }
         for (int k = 0; k < 2; ++k)
            if (br[k].Length())
               loci += TString::Format("%s{\"label\":\"%s%s\",\"pts\":[%s]}",
                                       loci.Length() > 1 ? "," : "", lab.Data(),
                                       k ? " (high-KE branch)" : "", br[k].Data());
         printf("  locus %-12s branches:%s%s   max lab angle %.0f deg\n", lab.Data(),
                br[0].Length() ? " low-KE" : "", br[1].Length() ? " high-KE" : "", aMax);
      }
   }
   loci += "]";

   TString gatesJson = "{";
   {
      TObjArray *toks = gatesCSV.Tokenize(",");
      bool first = true;
      for (int i = 0; i < toks->GetEntries(); ++i) {
         TString p = ((TObjString *)toks->At(i))->GetString().Strip(TString::kBoth);
         TString full = p.BeginsWith("/") ? p : here + "/" + p;
         std::string body = readAll(full.Data());
         if (body.empty()) { printf("  (gate not found, skipped: %s)\n", full.Data()); continue; }
         gatesJson += TString::Format("%s\"%s\":%s", first ? "" : ",",
                                      gSystem->BaseName(p.Data()), body.c_str());
         first = false;
         printf("  overlay gate: %s\n", full.Data());
      }
   }
   gatesJson += "}";

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
.panel{background:var(--pan);border:1px solid var(--line);border-radius:8px;padding:14px}
.panel h2{font-size:11px;text-transform:uppercase;letter-spacing:.12em;color:var(--dim);
margin:0 0 10px;font-weight:600}
label{display:block;margin:9px 0 3px;color:var(--dim);font-size:12px}
input,select,textarea{width:100%;background:#0d0f14;color:var(--fg);
border:1px solid var(--line);border-radius:5px;padding:6px 8px;font:inherit}
button{background:var(--acc);color:#12141a;border:0;border-radius:5px;padding:7px 11px;
font:inherit;font-weight:700;cursor:pointer;margin:3px 3px 0 0}
button.sec{background:#2a2f3a;color:var(--fg);font-weight:500}
canvas{border:1px solid var(--line);border-radius:8px;background:#0d0f14;display:block}
#cv{cursor:crosshair}
.stat{display:flex;gap:12px;flex-wrap:wrap;margin-top:10px}
.stat div{background:#0d0f14;border:1px solid var(--line);border-radius:6px;padding:8px 12px}
.stat b{display:block;font-size:17px;color:var(--acc)}
.stat span{font-size:11px;color:var(--dim);text-transform:uppercase;letter-spacing:.08em}
.hint{color:var(--dim);font-size:11.5px;margin-top:9px;line-height:1.6}
textarea{height:120px;font-size:11px}
.leg{display:flex;gap:12px;flex-wrap:wrap;font-size:11px;color:var(--dim);margin-top:6px}
.leg i{display:inline-block;width:14px;height:3px;vertical-align:middle;margin-right:4px}
</style></head><body>
<header><h1>AT-TPC PID gate editor</h1><div class="sub" id="sub"></div></header>
<div class="wrap">
 <div class="panel" style="flex:0 0 250px">
  <h2>display</h2>
  <label>track quality</label>
  <select id="q"><option value="0">all tracks</option>
  <option value="15" selected>nClusters &ge; 15</option><option value="25">nClusters &ge; 25</option></select>
  <label><input type="checkbox" id="logz" checked style="width:auto"> log colour scale</label>
  <label>overlay existing gate</label><select id="ov"></select>
  <h2 style="margin-top:16px">gate</h2>
  <button id="undo" class="sec">undo</button>
  <button id="clear" class="sec">clear</button>
  <button id="close">close</button>
  <button id="seed" class="sec">start from overlay</button>
  <label>name</label><input type="text" id="gname" value="triton_d2">
  <label>Z, A</label>
  <div style="display:flex;gap:6px"><input type="number" id="gz" value="1"><input type="number" id="ga" value="3"></div>
  <button id="dl">download JSON</button>
  <button id="exp" class="sec">show JSON</button>
  <div class="hint">Click to add vertices. The right panel shows B&rho; vs angle for the
  tracks currently inside &mdash; a good gate sits on its kinematic curve.</div>
 </div>
 <div>
  <canvas id="cv" width="620" height="480"></canvas>
  <div class="stat">
   <div><span>inside</span><b id="nin">0</b></div>
   <div><span>shown</span><b id="ntot">0</b></div>
   <div><span>fraction</span><b id="frac">0 %</b></div>
   <div><span>vertices</span><b id="nv">0</b></div>
  </div>
 </div>
 <div>
  <canvas id="cv2" width="620" height="480"></canvas>
  <div class="leg" id="leg"></div>
 </div>
 <div class="panel" style="flex:1 1 100%">
  <h2>json</h2><textarea id="out" readonly></textarea>
  <button id="copy" class="sec">copy</button>
  <div class="hint">Save into <code>pid/</code> and pass as the <code>pidGate</code> argument
  of the fit macro.</div>
 </div>
</div>
<script>
const XMAX=__XMAX__, YMAX=__YMAX__, NX=__NX__, NY=__NY__, TAG="__TAG__";
const GATES=__GATES__, LOCI=__LOCI__, EBEAM=__EBEAM__;
const RAW=Uint8Array.from(atob("__DATA__"), c=>c.charCodeAt(0));
const NT=RAW.length/4;
document.getElementById('sub').textContent =
  TAG+" · Spyral PID: √dE/dx vs Bρ · "+NT.toLocaleString()
  +" ungated pattern tracks · loci at E_beam = "+EBEAM+" MeV, Ex = 0";
const cv=document.getElementById('cv'),  cx=cv.getContext('2d');
const cv2=document.getElementById('cv2'),cx2=cv2.getContext('2d');
const L=64,R=14,T=14,B=46;
const PW=cv.width-L-R, PH=cv.height-T-B;
let poly=[], closed=false;
const ovSel=document.getElementById('ov');
ovSel.innerHTML='<option value="">none</option>'+
  Object.keys(GATES).map(k=>`<option value="${k}">${k}</option>`).join('');
const LCOL=['#e8933a','#6ec7ff','#7CFF9B','#ff7ba8','#c9a0ff'];
document.getElementById('leg').innerHTML = LOCI.map((l,i)=>
  `<span><i style="background:${LCOL[i%LCOL.length]}"></i>${l.label}</span>`).join('');

// decoded track values
const vx=new Float32Array(NT), vy=new Float32Array(NT), vp=new Float32Array(NT), vq=new Uint8Array(NT);
for(let i=0;i<NT;i++){ vx[i]=RAW[4*i]/255*XMAX; vy[i]=RAW[4*i+1]/255*YMAX;
                       vp[i]=RAW[4*i+2]/255*180; vq[i]=RAW[4*i+3]; }

function inPoly(x,y,P){let c=false;
  for(let i=0,j=P.length-1;i<P.length;j=i++){const a=P[i],b=P[j];
    if(((a[1]>y)!=(b[1]>y)) && (x<(b[0]-a[0])*(y-a[1])/(b[1]-a[1])+a[0])) c=!c;}
  return c;}
const ramp=f=>{f=Math.max(0,Math.min(1,f));
  const r=Math.round(255*Math.min(1,Math.max(0,-0.4+1.8*f-0.3*f*f)));
  const g=Math.round(255*Math.min(1,Math.max(0,0.05+1.15*f-0.35*f*f)));
  const b=Math.round(255*Math.min(1,Math.max(0,0.35+0.9*f-1.4*f*f)));
  return `rgb(${r},${g},${b})`;};

function axes(g,xlo,xhi,ylo,yhi,xlab,ylab,W,H){
  g.strokeStyle='#2a2f3a'; g.fillStyle='#8b93a7'; g.lineWidth=1;
  g.strokeRect(L,T,W,H); g.font='11px ui-monospace,monospace'; g.textAlign='center';
  for(let k=0;k<=8;k++){const v=xlo+(xhi-xlo)*k/8, px=L+W*k/8;
    g.beginPath();g.moveTo(px,T+H);g.lineTo(px,T+H+4);g.stroke();
    g.fillText(v.toFixed(v>=100?0:(xhi-xlo>10?0:2)),px,T+H+17);}
  g.fillText(xlab,L+W/2,T+H+36);
  g.textAlign='right';
  for(let k=0;k<=5;k++){const v=ylo+(yhi-ylo)*k/5, py=T+H-H*k/5;
    g.beginPath();g.moveTo(L,py);g.lineTo(L-4,py);g.stroke();
    g.fillText(v.toFixed(2),L-8,py+4);}
  g.save();g.translate(15,T+H/2);g.rotate(-Math.PI/2);g.textAlign='center';
  g.fillText(ylab,0,0);g.restore();}

function heat(g,pts,xlo,xhi,ylo,yhi,nbx,nby,logz){
  const h=new Int32Array(nbx*nby); let mx=0;
  for(const [px,py] of pts){
    const ix=Math.floor((px-xlo)/(xhi-xlo)*nbx), iy=Math.floor((py-ylo)/(yhi-ylo)*nby);
    if(ix<0||ix>=nbx||iy<0||iy>=nby) continue;
    const v=++h[iy*nbx+ix]; if(v>mx) mx=v;}
  const nz=[]; for(let i=0;i<h.length;i++) if(h[i]) nz.push(h[i]);
  nz.sort((a,b)=>a-b); const sat=nz.length?nz[Math.floor(nz.length*0.995)]:1;
  const bw=PW/nbx, bh=PH/nby;
  for(let iy=0;iy<nby;iy++)for(let ix=0;ix<nbx;ix++){
    const v=h[iy*nbx+ix]; if(!v) continue;
    const f=logz?Math.log10(1+v)/Math.log10(1+sat):v/sat;
    g.fillStyle=ramp(f);
    g.fillRect(L+ix*bw, T+PH-(iy+1)*bh, Math.ceil(bw), Math.ceil(bh));}
}

function draw(){
  const qmin=+document.getElementById('q').value;
  const logz=document.getElementById('logz').checked;
  const gx=v=>L+v/XMAX*PW, gy=v=>T+PH-v/YMAX*PH;

  // ---- left: the PID landscape ----
  cx.clearRect(0,0,cv.width,cv.height); cx.fillStyle='#0d0f14'; cx.fillRect(L,T,PW,PH);
  const land=[], sel=[];
  let nin=0, ntot=0;
  const usePoly = poly.length>2;
  for(let i=0;i<NT;i++){
    if(vq[i]<qmin) continue;
    land.push([vx[i],vy[i]]); ++ntot;
    if(usePoly && inPoly(vx[i],vy[i],poly)){ sel.push([vp[i],vy[i]]); ++nin; }
  }
  heat(cx,land,0,XMAX,0,YMAX,NX,NY,logz);
  axes(cx,0,XMAX,0,YMAX,'√dE/dx  [arb]','Bρ  [T·m]',PW,PH);
  const ov=ovSel.value;
  if(ov&&GATES[ov]){const V=GATES[ov].vertices;
    cx.strokeStyle='#6ec7ff';cx.lineWidth=1.5;cx.setLineDash([5,3]);cx.beginPath();
    V.forEach((p,i)=>i?cx.lineTo(gx(p[0]),gy(p[1])):cx.moveTo(gx(p[0]),gy(p[1])));
    cx.closePath();cx.stroke();cx.setLineDash([]);}
  if(poly.length){cx.strokeStyle='#7CFF9B';cx.lineWidth=2;cx.beginPath();
    poly.forEach((p,i)=>i?cx.lineTo(gx(p[0]),gy(p[1])):cx.moveTo(gx(p[0]),gy(p[1])));
    if(closed)cx.closePath(); cx.stroke();
    cx.fillStyle='#7CFF9B'; poly.forEach(p=>cx.fillRect(gx(p[0])-2.5,gy(p[1])-2.5,5,5));}

  document.getElementById('nin').textContent=nin.toLocaleString();
  document.getElementById('ntot').textContent=ntot.toLocaleString();
  document.getElementById('frac').textContent=ntot?(100*nin/ntot).toFixed(2)+' %':'0 %';
  document.getElementById('nv').textContent=poly.length;

  // ---- right: Brho vs polar angle for the SELECTED tracks ----
  cx2.clearRect(0,0,cv2.width,cv2.height); cx2.fillStyle='#0d0f14'; cx2.fillRect(L,T,PW,PH);
  if(sel.length) heat(cx2,sel,0,180,0,YMAX,128,128,logz);
  axes(cx2,0,180,0,YMAX,'θ_polar  [deg]','Bρ  [T·m]',PW,PH);
  LOCI.forEach((l,i)=>{
    cx2.strokeStyle=LCOL[i%LCOL.length]; cx2.lineWidth=2; cx2.beginPath();
    l.pts.forEach((p,k)=>{const X=L+p[0]/180*PW, Y=T+PH-p[1]/YMAX*PH;
      k?cx2.lineTo(X,Y):cx2.moveTo(X,Y);});
    cx2.stroke();});
  cx2.fillStyle='#8b93a7'; cx2.font='12px ui-monospace,monospace'; cx2.textAlign='left';
  cx2.fillText(usePoly?`selected: ${sel.length.toLocaleString()} tracks`
                      :'draw a gate to populate this panel', L+8, T+16);
}

const px2x=p=>(p-L)/PW*XMAX, py2y=p=>(T+PH-p)/PH*YMAX;
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
['q','logz','ov'].forEach(id=>document.getElementById(id).addEventListener('input',draw));
function json(){
  const o={name:document.getElementById('gname').value,xaxis:"sqrtdedx",yaxis:"brho",vertices:poly};
  const z=+document.getElementById('gz').value, a=+document.getElementById('ga').value;
  if(z) o.Z=z; if(a) o.A=a;
  return JSON.stringify(o,null,1);}
document.getElementById('exp').onclick=()=>{document.getElementById('out').value=json();};
document.getElementById('copy').onclick=()=>{const t=document.getElementById('out');
  t.value=t.value||json(); t.select(); document.execCommand('copy');};
document.getElementById('dl').onclick=()=>{
  if(poly.length<3){alert('draw at least 3 vertices first');return;}
  document.getElementById('out').value=json();
  const b=new Blob([json()],{type:'application/json'});
  const u=URL.createObjectURL(b), a=document.createElement('a');
  a.href=u; a.download=document.getElementById('gname').value+'.json';
  document.body.appendChild(a); a.click(); document.body.removeChild(a);
  setTimeout(()=>URL.revokeObjectURL(u),4000);};
draw();
</script></body></html>)HTML");
   fclose(o);

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
   rep("__EBEAM__", TString::Format("%g", Ebeam).Data());
   rep("__GATES__", gatesJson.Data());
   rep("__LOCI__", loci.Data());
   rep("__DATA__", b64(rec));
   o = fopen(outHtml.Data(), "w");
   fwrite(html.data(), 1, html.size(), o);
   fclose(o);
   printf("wrote %s  (%.1f MB)\n", outHtml.Data(), html.size()/1e6);
}
