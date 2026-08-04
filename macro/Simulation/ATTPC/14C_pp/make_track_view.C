/// @file make_track_view.C
/// @brief Dump reconstructed events (hits + PRA track assignment) into a self-contained
/// HTML point-cloud viewer: three projections, event navigation, colour by track.
///
/// ROOT's Eve viewer crashes under WSL (no OpenGL), so this goes to the browser instead --
/// the same route the excitation explorer already uses.
///
///   root -b -q 'make_track_view.C("./data/sim_disp.root","/home/yassid/c14_tracks.html",60)'

void make_track_view(TString inFile = "./data/sim_disp.root", TString outHtml = "/home/yassid/c14_tracks.html",
                     Int_t maxEvents = 60, Int_t minHits = 8)
{
   gSystem->Load("libAtReconstruction.so");
   TFile *f = TFile::Open(inFile);
   if (!f || f->IsZombie()) {
      printf("\033[1;31mcannot open %s\033[0m\n", inFile.Data());
      return;
   }
   auto *t = (TTree *)f->Get("cbmsim");
   TClonesArray *ev = nullptr, *pe = nullptr;
   t->SetBranchAddress("AtEventH", &ev);
   t->SetBranchAddress("AtPatternEvent", &pe);

   FILE *o = fopen(outHtml.Data(), "w");
   if (!o) {
      printf("cannot write %s\n", outHtml.Data());
      return;
   }

   // ---- collect events -----------------------------------------------------
   std::string js = "[";
   int nOut = 0;
   for (Long64_t i = 0; i < t->GetEntries() && nOut < maxEvents; ++i) {
      t->GetEntry(i);
      if (!ev || !ev->GetEntriesFast())
         continue;
      auto *e = (AtEvent *)ev->At(0);
      if (!e || e->GetNumHits() < minHits)
         continue;

      // hit -> track id from the pattern event (-1 = not on any track)
      std::map<int, int> hitTrack;
      int nTracks = 0;
      if (pe && pe->GetEntriesFast()) {
         auto *p = (AtPatternEvent *)pe->At(0);
         if (p) {
            int tid = 0;
            for (auto &tr : p->GetTrackCand()) {
               for (auto &h : tr.GetHitArray()) // HitVector holds pointers
                  if (h)
                     hitTrack[h->GetHitID()] = tid;
               ++tid;
            }
            nTracks = tid;
         }
      }
      if (nTracks == 0)
         continue;

      std::string hx = "[", hy = "[", hz = "[", hq = "[", ht = "[";
      int n = 0;
      for (int h = 0; h < e->GetNumHits(); ++h) {
         auto &hit = e->GetHit(h);
         auto pos = hit.GetPosition();
         int tid = hitTrack.count(hit.GetHitID()) ? hitTrack[hit.GetHitID()] : -1;
         if (n) {
            hx += ",";
            hy += ",";
            hz += ",";
            hq += ",";
            ht += ",";
         }
         hx += TString::Format("%.1f", pos.X()).Data();
         hy += TString::Format("%.1f", pos.Y()).Data();
         hz += TString::Format("%.1f", pos.Z()).Data();
         hq += TString::Format("%.0f", hit.GetCharge()).Data();
         ht += TString::Format("%d", tid).Data();
         ++n;
      }
      hx += "]";
      hy += "]";
      hz += "]";
      hq += "]";
      ht += "]";
      if (nOut)
         js += ",";
      js += TString::Format("{\"ev\":%lld,\"nt\":%d,\"x\":%s,\"y\":%s,\"z\":%s,\"q\":%s,\"t\":%s}", i, nTracks,
                            hx.c_str(), hy.c_str(), hz.c_str(), hq.c_str(), ht.c_str())
               .Data();
      ++nOut;
   }
   js += "]";

   // fputs, NOT fprintf: the embedded JS contains '%' (modulo) which fprintf would eat
   // as format specifiers.
   fputs(R"HTML(<meta charset="utf-8"><title>14C sim - tracks</title>
<style>
:root{--bg:#0e1116;--pan:#161b22;--ln:#30363d;--ink:#e6edf3;--ink2:#8b949e;--acc:#f0883e}
@media(prefers-color-scheme:light){:root{--bg:#fff;--pan:#f6f8fa;--ln:#d0d7de;--ink:#1f2328;--ink2:#656d76}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);
font:13px/1.45 ui-monospace,SFMono-Regular,Menlo,monospace}
header{padding:10px 16px;border-bottom:1px solid var(--ln);display:flex;gap:14px;align-items:center;flex-wrap:wrap}
h1{font-size:15px;margin:0;font-weight:600}
button{background:var(--pan);color:var(--ink);border:1px solid var(--ln);border-radius:6px;
padding:5px 11px;cursor:pointer;font:inherit}button:hover{border-color:var(--acc)}
.wrap{display:grid;grid-template-columns:repeat(auto-fit,minmax(340px,1fr));gap:12px;padding:12px}
.card{background:var(--pan);border:1px solid var(--ln);border-radius:9px;padding:9px}
.card h2{font-size:12px;margin:0 0 6px;font-weight:600;color:var(--ink2)}
canvas{width:100%;height:300px;display:block}
#info{color:var(--ink2)}
</style>
<header>
  <h1>14C(p,p') simulation &mdash; reconstructed hits</h1>
  <button id="prev">&larr; prev</button><button id="next">next &rarr;</button>
  <span id="info"></span>
  <label style="color:var(--ink2)"><input type="checkbox" id="onlyTrk"> only hits on tracks</label>
</header>
<div class="wrap">
  <div class="card"><h2>x &ndash; z  (beam along z)</h2><canvas id="cxz"></canvas></div>
  <div class="card"><h2>y &ndash; z</h2><canvas id="cyz"></canvas></div>
  <div class="card"><h2>x &ndash; y  (pad plane)</h2><canvas id="cxy"></canvas></div>
</div>
<script>
const EV = )HTML",
         o);
   fputs(js.c_str(), o);
   fputs(R"HTML(;
let idx = 0;
const $ = i => document.getElementById(i);
const PAL = ['#f0883e','#3fb950','#58a6ff','#d29922','#db61a2','#a371f7','#2ea043','#f85149'];
const css = v => getComputedStyle(document.documentElement).getPropertyValue(v).trim();

function draw(cv, xs, ys, ts, xl, yl){
  const dpr = devicePixelRatio||1, w = cv.clientWidth, h = cv.clientHeight;
  cv.width = w*dpr; cv.height = h*dpr;
  const g = cv.getContext('2d'); g.setTransform(dpr,0,0,dpr,0,0);
  g.clearRect(0,0,w,h);
  const P = {l:44,r:10,t:8,b:26}, pw = w-P.l-P.r, ph = h-P.t-P.b;
  if(!xs.length){return;}
  let x0=Math.min(...xs), x1=Math.max(...xs), y0=Math.min(...ys), y1=Math.max(...ys);
  const padx=(x1-x0)*0.06+1, pady=(y1-y0)*0.06+1;
  x0-=padx; x1+=padx; y0-=pady; y1+=pady;
  g.strokeStyle=css('--ln'); g.lineWidth=1; g.strokeRect(P.l,P.t,pw,ph);
  g.fillStyle=css('--ink2'); g.font='10px ui-monospace';
  g.textAlign='right'; g.textBaseline='middle';
  for(let k=0;k<=4;k++){ const v=y0+(y1-y0)*k/4, yy=P.t+ph-ph*k/4;
    g.fillText(v.toFixed(0),P.l-5,yy); }
  g.textAlign='center'; g.textBaseline='top';
  for(let k=0;k<=4;k++){ const v=x0+(x1-x0)*k/4, xx=P.l+pw*k/4;
    g.fillText(v.toFixed(0),xx,P.t+ph+5); }
  g.textAlign='right'; g.fillText(xl,P.l+pw,P.t+ph+13);
  g.save(); g.translate(11,P.t); g.rotate(-Math.PI/2); g.textAlign='right'; g.textBaseline='top';
  g.fillText(yl,0,0); g.restore();
  for(let i=0;i<xs.length;i++){
    const px=P.l+(xs[i]-x0)/(x1-x0)*pw, py=P.t+ph-(ys[i]-y0)/(y1-y0)*ph;
    g.fillStyle = ts[i]<0 ? 'rgba(139,148,158,.35)' : PAL[ts[i]%PAL.length];
    g.beginPath(); g.arc(px,py,ts[i]<0?1.4:2.4,0,6.283); g.fill();
  }
}
function render(){
  if(!EV.length){ $('info').textContent='no events with tracks'; return; }
  const e = EV[idx], only = $('onlyTrk').checked;
  const keep = [];
  for(let i=0;i<e.x.length;i++) if(!only || e.t[i]>=0) keep.push(i);
  const X=keep.map(i=>e.x[i]), Y=keep.map(i=>e.y[i]), Z=keep.map(i=>e.z[i]), T=keep.map(i=>e.t[i]);
  draw($('cxz'),Z,X,T,'z [mm]','x [mm]');
  draw($('cyz'),Z,Y,T,'z [mm]','y [mm]');
  draw($('cxy'),X,Y,T,'x [mm]','y [mm]');
  const onTrk = e.t.filter(v=>v>=0).length;
  $('info').textContent = `event ${e.ev}  (${idx+1}/${EV.length})   hits ${e.x.length}   on tracks ${onTrk}   tracks ${e.nt}`;
}
$('prev').onclick=()=>{idx=(idx-1+EV.length)%EV.length;render()};
$('next').onclick=()=>{idx=(idx+1)%EV.length;render()};
$('onlyTrk').onchange=render;
addEventListener('keydown',e=>{if(e.key==='ArrowLeft')$('prev').click();if(e.key==='ArrowRight')$('next').click()});
addEventListener('resize',render);
matchMedia('(prefers-color-scheme:dark)').addEventListener('change',render);
render();
</script>
)HTML",
         o);
   fclose(o);
   f->Close();
   printf("\nwrote %s   (%d events)\n", outHtml.Data(), nOut);
}
