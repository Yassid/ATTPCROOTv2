/// @file tsync_check.C
/// @brief Offset-independent alignment check between pad stream (AtEventH.fTimestamp) and
///        FRIB stream, using timestamp DELTAS (pad[i+1]-pad[i] vs frib[i+1]-frib[i]).
///        Reports how many events are aligned vs the first real desync (a dropped event),
///        which entry-index IC-matching would get wrong.
///   root -b -q 'tsync_check.C("run_0143")'
void tsync_check(TString run, TString recoDir = "/home/yassid/a2091_C15_reco/",
                 TString fribDir = "/home/yassid/a2091_C15_reco_slim/")
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
   // The two DAQ clocks do NOT always tick at the same rate in a2091: in the early 15C runs
   // pad ts1 shares the FRIB clock (scale 1), while from ~run_0138 on both pad timestamps run
   // 100x finer than the FRIB one. So scan the scale as well as the index, else an aligned run
   // reads as 0% agreement. Agreement is on dPad ~ scale*dFrib within tol (in FRIB ticks).
   auto agree = [&](std::vector<Long64_t> &p, double sc, double tol) {
      long ok = 0, tot = 0;
      for (size_t i = 1; i < N; i++) {
         if (p[i] < 0 || p[i - 1] < 0 || fr[i] < 0 || fr[i - 1] < 0) continue;
         ++tot;
         if (std::fabs((double)(p[i] - p[i - 1]) / sc - (double)(fr[i] - fr[i - 1])) <= tol) ++ok;
      }
      return tot ? 100.0 * ok / tot : 0.0;
   };
   const double scales[] = {1.0, 100.0};
   int bIdx = 1; double bSc = 1.0, bAg = -1;
   for (int idx = 0; idx < 2; ++idx)
      for (double sc : scales) {
         double a = agree(idx ? p1 : p0, sc, 2.0);
         if (a > bAg) { bAg = a; bIdx = idx; bSc = sc; }
      }
   auto &pad = (bIdx ? p1 : p0);
   printf("%s  padN=%zu fribN=%zu (diff %+ld)  idx%d scale%.0f  delta-agree tol0/2/5/20/100 = "
          "%.1f/%.1f/%.1f/%.1f/%.1f%%\n",
          run.Data(), p0.size(), fr.size(), (long)fr.size() - (long)p0.size(), bIdx, bSc, agree(pad, bSc, 0),
          agree(pad, bSc, 2), agree(pad, bSc, 5), agree(pad, bSc, 20), agree(pad, bSc, 100));
}
