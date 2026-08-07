/// @file make_digi_view.C
/// @brief Browser viewer for DIGITIZED point clouds: three projections, event navigation,
///        colour by PRA track or by charge, and an A/B switch to flip between two files
///        (e.g. the simulation and the a1954 data) without leaving the page.
///
/// Works on either kind of reco file:
///   * sim reco      -- has AtEventH, so ALL digitized hits are shown and the ones that no
///                      PRA track claimed are drawn in grey. That is the view that shows
///                      fragmentation directly.
///   * data gated    -- ~/a1954_C14_fit_300torr/in/ carries AtPatternEvent only, so only
///                      on-track hits exist. The page says which mode each dataset is in.
///
/// ROOT's Eve viewer cannot run under WSL (no OpenGL), hence the browser route.
///
///   root -b -q 'make_digi_view.C("fileA.root","SIM","fileB.root","DATA","/home/yassid/c14_digi.html",40)'

static std::string dumpFile(TString inFile, TString label, Int_t maxEvents, Int_t minHits, bool &ok)
{
   ok = false;
   TFile *f = TFile::Open(inFile);
   if (!f || f->IsZombie()) {
      printf("\033[1;31mcannot open %s\033[0m\n", inFile.Data());
      return "";
   }
   auto *t = (TTree *)f->Get("cbmsim");
   if (!t) {
      printf("\033[1;31mno cbmsim in %s\033[0m\n", inFile.Data());
      return "";
   }
   const bool hasEvt = t->GetBranch("AtEventH") != nullptr;
   const bool hasPat = t->GetBranch("AtPatternEvent") != nullptr;
   if (!hasPat && !hasEvt) {
      printf("\033[1;31m%s has neither AtEventH nor AtPatternEvent\033[0m\n", inFile.Data());
      return "";
   }
   TClonesArray *ev = nullptr, *pe = nullptr;
   if (hasEvt)
      t->SetBranchAddress("AtEventH", &ev);
   if (hasPat)
      t->SetBranchAddress("AtPatternEvent", &pe);

   std::string js = "[";
   int nOut = 0;
   for (Long64_t i = 0; i < t->GetEntries() && nOut < maxEvents; ++i) {
      t->GetEntry(i);

      // hit id -> track index, from the pattern event (-1 = claimed by no track)
      std::map<int, int> hitTrack;
      int nTracks = 0;
      AtPatternEvent *p = (hasPat && pe && pe->GetEntriesFast()) ? (AtPatternEvent *)pe->At(0) : nullptr;
      if (p) {
         int tid = 0;
         for (auto &tr : p->GetTrackCand()) {
            for (auto &h : tr.GetHitArray()) // HitVector holds POINTERS
               if (h)
                  hitTrack[h->GetHitID()] = tid;
            ++tid;
         }
         nTracks = tid;
      }

      std::string hx = "[", hy = "[", hz = "[", hq = "[", ht = "[", hm = "[";
      int n = 0;
      // mc: MC truth trackID of the hit (-1 unlabelled), so the page can colour by TRUTH and
      // the clustering can be judged by eye against it. Present only when MC info was saved.
      auto push = [&](double x, double y, double z, double q, int tid, int mcid) {
         if (n) { hx += ","; hy += ","; hz += ","; hq += ","; ht += ","; hm += ","; }
         hx += TString::Format("%.1f", x).Data();
         hy += TString::Format("%.1f", y).Data();
         hz += TString::Format("%.1f", z).Data();
         hq += TString::Format("%.0f", q).Data();
         ht += TString::Format("%d", tid).Data();
         hm += TString::Format("%d", mcid).Data();
         ++n;
      };

      std::set<int> mcSeen;
      if (hasEvt && ev && ev->GetEntriesFast()) {
         auto *e = (AtEvent *)ev->At(0);
         if (e)
            for (int h = 0; h < e->GetNumHits(); ++h) {
               auto &hit = e->GetHit(h);
               auto pos = hit.GetPosition();
               const auto &mc = hit.GetMCSimPointArray();
               int mcid = mc.empty() ? -1 : mc[0].trackID;
               if (mcid >= 0)
                  mcSeen.insert(mcid);
               push(pos.X(), pos.Y(), pos.Z(), hit.GetCharge(),
                    hitTrack.count(hit.GetHitID()) ? hitTrack[hit.GetHitID()] : -1, mcid);
            }
      } else if (p) { // data: on-track hits only, no MC truth
         int tid = 0;
         for (auto &tr : p->GetTrackCand()) {
            for (auto &h : tr.GetHitArray())
               if (h) {
                  auto pos = h->GetPosition();
                  push(pos.X(), pos.Y(), pos.Z(), h->GetCharge(), tid, -1);
               }
            ++tid;
         }
      }
      if (n < minHits)
         continue;

      hx += "]"; hy += "]"; hz += "]"; hq += "]"; ht += "]"; hm += "]";
      if (nOut)
         js += ",";
      js += TString::Format("{\"ev\":%lld,\"nt\":%d,\"nm\":%d,\"x\":%s,\"y\":%s,\"z\":%s,\"q\":%s,\"t\":%s,\"m\":%s}",
                            i, nTracks, (int)mcSeen.size(), hx.c_str(), hy.c_str(), hz.c_str(), hq.c_str(),
                            ht.c_str(), hm.c_str())
               .Data();
      ++nOut;
   }
   js += "]";
   f->Close();
   printf("  %-6s %-60s  %d events, %s\n", label.Data(), gSystem->BaseName(inFile), nOut,
          hasEvt ? "ALL digitized hits (off-track shown grey)" : "on-track hits only");
   ok = (nOut > 0);
   return js;
}

void make_digi_view(TString fileA = "./data/sim_reco.root", TString labelA = "SIM",
                    TString fileB = "/home/yassid/a1954_C14_fit_300torr/in/run_0056_reco.root",
                    TString labelB = "DATA", TString outHtml = "/home/yassid/c14_digi.html", Int_t maxEvents = 40,
                    Int_t minHits = 8, TString fileC = "", TString labelC = "C")
{
   gSystem->Load("libAtReconstruction.so");
   printf("\ncollecting:\n");
   bool okA = false, okB = false, okC = false;
   std::string jsA = dumpFile(fileA, labelA, maxEvents, minHits, okA);
   std::string jsB = fileB.Length() ? dumpFile(fileB, labelB, maxEvents, minHits, okB) : "[]";
   std::string jsC = fileC.Length() ? dumpFile(fileC, labelC, maxEvents, minHits, okC) : "[]";
   if (!okA && !okB && !okC) {
      printf("\033[1;31mnothing to show\033[0m\n");
      return;
   }
   if (!okB)
      jsB = "[]";
   if (!okC)
      jsC = "[]";

   FILE *o = fopen(outHtml.Data(), "w");
   if (!o) {
      printf("cannot write %s\n", outHtml.Data());
      return;
   }

   fputs("<!doctype html><html><head><meta charset='utf-8'><title>14C digitized tracks</title>\n", o);
   fputs("<style>\n"
         ":root{--bg:#0e1116;--fg:#e6edf3;--dim:#8b949e;--card:#161b22;--line:#30363d;--accent:#f0883e}\n"
         "@media(prefers-color-scheme:light){:root{--bg:#fff;--fg:#1f2328;--dim:#656d76;--card:#f6f8fa;"
         "--line:#d0d7de;--accent:#bc4c00}}\n"
         "*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--fg);"
         "font:13px ui-monospace,SFMono-Regular,Menlo,monospace}\n"
         "header{padding:14px 18px;border-bottom:1px solid var(--line);display:flex;gap:18px;align-items:center;"
         "flex-wrap:wrap}\n"
         "h1{font-size:15px;margin:0;font-weight:600}\n"
         ".sub{color:var(--dim);font-size:11px}\n"
         "button{background:var(--card);color:var(--fg);border:1px solid var(--line);border-radius:6px;"
         "padding:6px 12px;cursor:pointer;font:inherit}\n"
         "button:hover{border-color:var(--accent)}button.on{background:var(--accent);color:#fff;"
         "border-color:var(--accent)}\n"
         "input{background:var(--card);color:var(--fg);border:1px solid var(--line);border-radius:6px;"
         "padding:5px 8px;width:70px;font:inherit}\n"
         ".wrap{display:grid;grid-template-columns:repeat(auto-fit,minmax(420px,1fr));gap:14px;padding:14px}\n"
         ".card{background:var(--card);border:1px solid var(--line);border-radius:10px;padding:10px}\n"
         ".ttl{font-size:11px;color:var(--dim);margin-bottom:6px}\n"
         "canvas{width:100%;height:auto;display:block}\n"
         ".stat{color:var(--dim)}.stat b{color:var(--fg)}\n"
         "</style></head><body>\n",
         o);
   fputs("<header><h1>14C(p,p&#39;) &middot; digitized tracks</h1>\n", o);
   fprintf(o, "<div><button id='dsA' class='on'>%s</button><button id='dsB'>%s</button>%s</div>\n", labelA.Data(),
           labelB.Data(),
           okC ? TString::Format("<button id='dsC'>%s</button>", labelC.Data()).Data() : "");
   fputs("<div><button id='prev'>&larr; prev</button> <input id='jump' type='number' value='0'> "
         "<button id='next'>next &rarr;</button></div>\n"
         "<div><button id='cTrk' class='on'>colour: cluster</button><button id='cMC'>colour: truth</button>"
         "<button id='cQ'>colour: charge</button></div>\n"
         "<div class='sub'>point size <input id='psz' type='number' value='2.6' step='0.4' min='0.6'></div>\n"
         "<div class='stat' id='stat'></div></header>\n"
         "<div class='wrap'>\n"
         "<div class='card'><div class='ttl'>pad plane &nbsp; x &ndash; y</div><canvas id='cxy'></canvas></div>\n"
         "<div class='card'><div class='ttl'>drift &nbsp; z &ndash; x</div><canvas id='czx'></canvas></div>\n"
         "<div class='card'><div class='ttl'>drift &nbsp; z &ndash; y</div><canvas id='czy'></canvas></div>\n"
         "</div>\n<script>\n",
         o);
   fprintf(o, "const DS={A:{name:%s,ev:%s},B:{name:%s,ev:%s},C:{name:%s,ev:%s}};\n", ("\"" + labelA + "\"").Data(),
           jsA.c_str(), ("\"" + labelB + "\"").Data(), jsB.c_str(), ("\"" + labelC + "\"").Data(), jsC.c_str());
   fputs(
      "let cur='A', idx=0, mode='trk';\n"
      "const PAL=['#f0883e','#58a6ff','#3fb950','#db61a2','#e3b341','#a371f7','#39c5cf','#ff7b72'];\n"
      "function evs(){return DS[cur].ev}\n"
      "function qcol(v,lo,hi){const t=Math.max(0,Math.min(1,(v-lo)/(hi-lo||1)));\n"
      "  const r=Math.round(30+225*t), g=Math.round(60+120*(1-Math.abs(t-.5)*2)), b=Math.round(220-190*t);\n"
      "  return `rgb(${r},${g},${b})`}\n"
      "function draw(cv,X,Y,e,xlab,ylab){\n"
      "  const dpr=window.devicePixelRatio||1, W=cv.clientWidth||440, H=Math.round(W*0.72);\n"
      "  cv.width=W*dpr; cv.height=H*dpr; const c=cv.getContext('2d'); c.setTransform(dpr,0,0,dpr,0,0);\n"
      "  const cs=getComputedStyle(document.body); c.clearRect(0,0,W,H);\n"
      "  const m=34, iw=W-m-10, ih=H-m-10;\n"
      "  let x0=Math.min(...X),x1=Math.max(...X),y0=Math.min(...Y),y1=Math.max(...Y);\n"
      "  const px=(x1-x0)*0.06+1, py=(y1-y0)*0.06+1; x0-=px;x1+=px;y0-=py;y1+=py;\n"
      "  const sx=v=>m+(v-x0)/(x1-x0)*iw, sy=v=>H-m-(v-y0)/(y1-y0)*ih;\n"
      "  c.strokeStyle=cs.getPropertyValue('--line'); c.lineWidth=1;\n"
      "  c.strokeRect(m,H-m-ih,iw,ih);\n"
      "  c.fillStyle=cs.getPropertyValue('--dim'); c.font='10px ui-monospace,monospace';\n"
      "  c.fillText(xlab,W/2-14,H-8); c.save(); c.translate(11,H/2+14); c.rotate(-Math.PI/2);\n"
      "  c.fillText(ylab,0,0); c.restore();\n"
      "  c.fillText(x0.toFixed(0),m,H-m+13); c.fillText(x1.toFixed(0),m+iw-22,H-m+13);\n"
      "  c.fillText(y0.toFixed(0),4,H-m); c.fillText(y1.toFixed(0),4,H-m-ih+9);\n"
      "  const q=e.q, lo=Math.min(...q), hi=Math.max(...q);\n"
      "  const r=parseFloat(document.getElementById('psz').value)||2.6;\n"
      "  for(let i=0;i<X.length;i++){\n"
      "    const t = mode==='mc' ? e.m[i] : e.t[i];\n"
      "    c.fillStyle = mode==='q' ? qcol(q[i],lo,hi) : (t<0? '#6e7681' : PAL[t%PAL.length]);\n"
      "    c.globalAlpha = (mode!=='q'&&t<0)?0.45:0.9;\n"
      "    c.beginPath(); c.arc(sx(X[i]),sy(Y[i]),r,0,6.283); c.fill();}\n"
      "  c.globalAlpha=1;\n"
      "}\n"
      "function render(){\n"
      "  const E=evs(); const st=document.getElementById('stat');\n"
      "  if(!E.length){st.innerHTML='<b>no events in this dataset</b>';\n"
      "    ['cxy','czx','czy'].forEach(id=>{const cv=document.getElementById(id);\n"
      "      const c=cv.getContext('2d'); c.clearRect(0,0,cv.width,cv.height);}); return;}\n"
      "  idx=(idx%E.length+E.length)%E.length; const e=E[idx];\n"
      "  document.getElementById('jump').value=idx;\n"
      "  const off=e.t.filter(v=>v<0).length;\n"
      "  st.innerHTML=`<b>${DS[cur].name}</b> &nbsp; event ${idx+1}/${E.length} (entry ${e.ev}) &nbsp; `+\n"
      "    `hits <b>${e.x.length}</b> &nbsp; clusters <b>${e.nt}</b>`+(e.nm?` &nbsp; true tracks <b>${e.nm}</b>`:'')+\n"
      "    (off?` &nbsp; off-track <b>${off}</b>`:'');\n"
      "  draw(document.getElementById('cxy'),e.x,e.y,e,'x [mm]','y [mm]');\n"
      "  draw(document.getElementById('czx'),e.z,e.x,e,'z [mm]','x [mm]');\n"
      "  draw(document.getElementById('czy'),e.z,e.y,e,'z [mm]','y [mm]');\n"
      "}\n"
      "document.getElementById('prev').onclick=()=>{idx--;render()};\n"
      "document.getElementById('next').onclick=()=>{idx++;render()};\n"
      "document.getElementById('jump').oninput=ev=>{idx=parseInt(ev.target.value)||0;render()};\n"
      "document.getElementById('psz').oninput=render;\n"
      "function setDS(d){cur=d;idx=0;['A','B','C'].forEach(k=>{const b=document.getElementById('ds'+k);\n"
      "  if(b)b.classList.toggle('on',k===d)});render()}\n"
      "document.getElementById('dsA').onclick=()=>setDS('A');\n"
      "document.getElementById('dsB').onclick=()=>setDS('B');\n"
      "if(document.getElementById('dsC'))document.getElementById('dsC').onclick=()=>setDS('C');\n"
      "function setMode(m){mode=m;['trk','mc','q'].forEach(k=>{const b=document.getElementById(\n"
      "  k==='trk'?'cTrk':k==='mc'?'cMC':'cQ'); if(b)b.classList.toggle('on',k===m)});render()}\n"
      "document.getElementById('cTrk').onclick=()=>setMode('trk');\n"
      "document.getElementById('cMC').onclick=()=>setMode('mc');\n"
      "document.getElementById('cQ').onclick=()=>setMode('q');\n"
      "addEventListener('keydown',e=>{if(e.key==='c')setMode(mode==='trk'?'mc':'trk')});\n"
      "addEventListener('keydown',e=>{if(e.key==='ArrowLeft'){idx--;render()}if(e.key==='ArrowRight'){idx++;render()}});\n"
      "addEventListener('resize',render); render();\n"
      "</script></body></html>\n",
      o);
   fclose(o);

   FileStat_t st;
   gSystem->GetPathInfo(outHtml, st);
   printf("\nwrote %s  (%.1f MB)\n\n", outHtml.Data(), st.fSize / 1048576.);
}
