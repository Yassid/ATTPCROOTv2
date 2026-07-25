/// @file tsync_check.C
/// @brief Offset-independent alignment check between pad stream (AtEventH.fTimestamp) and
///        FRIB stream, using timestamp DELTAS (pad[i+1]-pad[i] vs frib[i+1]-frib[i]).
///        Reports how many events are aligned vs the first real desync (a dropped event),
///        which entry-index IC-matching would get wrong.
///   root -b -q 'tsync_check.C("run_0143")'
void tsync_check(TString run, TString recoDir = "/mnt/f/a1954_Be12_reco_hdb/",
                 TString fribDir = "/home/yassid/a1954_Be12_reco_hdb_slim/")
{
   gSystem->Load("libAtReconstruction.so");
   // pad timestamps (both indices), cheap: only timestamp leaf active
   TFile *fR = TFile::Open(recoDir + run + "_reco.root");
   TTree *tR = (TTree *)fR->Get("cbmsim");
   tR->SetBranchStatus("*", 0); tR->SetBranchStatus("AtEventH.fTimestamp*", 1);
   TClonesArray *ev = nullptr; tR->SetBranchAddress("AtEventH", &ev);
   std::vector<Long64_t> p0, p1;
   for (Long64_t i = 0; i < tR->GetEntries(); i++) { tR->GetEntry(i);
      if (ev->GetEntries()==0){p0.push_back(-1);p1.push_back(-1);continue;}
      auto *e=(AtEvent*)ev->At(0); p0.push_back((Long64_t)e->GetTimestamp(0));
      p1.push_back(e->GetTimestamps().size()>1?(Long64_t)e->GetTimestamp(1):-1); }
   fR->Close();
   TFile *fF = TFile::Open(fribDir + run + "_FRIB.root");
   TTree *tF = (TTree *)fF->Get("cbmsim");
   TClonesArray *fa = nullptr; tF->SetBranchAddress("AtRawEvent", &fa);
   std::vector<Long64_t> fr;
   for (Long64_t i = 0; i < tF->GetEntries(); i++) { tF->GetEntry(i);
      if (fa->GetEntries()==0){fr.push_back(-1);continue;} auto *r=(AtRawEvent*)fa->At(0); fr.push_back((Long64_t)r->GetTimestamp(0)); }
   fF->Close();

   size_t N = std::min({p0.size(), p1.size(), fr.size()});
   // pick pad index whose DELTAS best match frib deltas
   auto deltaAgree = [&](std::vector<Long64_t>&p){ long ok=0,tot=0; for(size_t i=1;i<N;i++){ if(p[i]<0||p[i-1]<0||fr[i]<0||fr[i-1]<0)continue; ++tot; if((p[i]-p[i-1])==(fr[i]-fr[i-1]))++ok; } return std::make_pair(ok,tot); };
   auto a0=deltaAgree(p0), a1=deltaAgree(p1);
   int idx = (a1.first>=a0.first)?1:0; auto &pad=(idx?p1:p0);
   // tolerance sweep on |dPad - dFrib|
   auto agreeTol=[&](long tol){ long ok=0,tot=0; for(size_t i=1;i<N;i++){ if(pad[i]<0||pad[i-1]<0||fr[i]<0||fr[i-1]<0)continue; ++tot; if(std::llabs((pad[i]-pad[i-1])-(fr[i]-fr[i-1]))<=tol)++ok; } return tot?100.0*ok/tot:0; };
   printf("%s  padN=%zu fribN=%zu (diff %+ld)  idx%d  delta-agree tol0/2/5/20/100 = %.1f/%.1f/%.1f/%.1f/%.1f%%\n",
          run.Data(), p0.size(), fr.size(), (long)fr.size()-(long)p0.size(), idx,
          agreeTol(0),agreeTol(2),agreeTol(5),agreeTol(20),agreeTol(100));
}
