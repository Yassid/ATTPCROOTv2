/// @file testMerge_a1975.C
/// @brief Validate AtGenfitter continuity merging: fit the same events with merging
/// OFF and ON, report input-track / fitted-track counts and how often merging fired.
///   root -l -b -q 'testMerge_a1975.C("run_0106", 150, "/mnt/f/a1975/reco/")'

// mode: 0=off, 1=default criteria, 2=ultra-loose (force any in-event merge)
void runPass(const char *tag, int mode, TString recoFile, Long64_t nEv, int &nInTracks, int &nFitted, int &nMergeEvents,
             int &nTracksAfter, int &nMultiTrackEvents)
{
   bool merge = (mode != 0);
   auto *f = TFile::Open(recoFile);
   auto *t = (TTree *)f->Get("cbmsim");
   auto *peArr = new TClonesArray("AtPatternEvent");
   t->SetBranchAddress("AtPatternEvent", &peArr);

   // proton genfit hypothesis (same as fitBoth), optional continuity merge
   auto fitter = std::make_unique<EventFit::AtGenfitter>(-2.85, 2212, 938.27208816 / 931.49410242, 1, "", kTRUE, 2, 5);
   fitter->SetZPadPlane(1000.0);
   fitter->SetMeasSigma(4.0);
   if (mode == 1)
      fitter->SetMergeContinuity(kTRUE, 50.0, 0.3, 30.0);
   else if (mode == 2)
      fitter->SetMergeContinuity(kTRUE, 1e9, 1e9, 1e9); // force: merge every pair in an event
   fitter->Init();

   nInTracks = nFitted = nMergeEvents = nTracksAfter = nMultiTrackEvents = 0;
   Long64_t n = std::min(nEv, t->GetEntries());
   for (Long64_t i = 0; i < n; ++i) {
      t->GetEntry(i);
      auto *pe = (AtPatternEvent *)peArr->At(0);
      if (!pe)
         continue;
      int nin = pe->GetTrackCand().size();
      nInTracks += nin;
      if (nin >= 2)
         nMultiTrackEvents++;
      auto te = std::make_unique<AtTrackingEvent>();
      fitter->FitEvent(te.get(), pe);
      int nAfter = te->GetTrackArray().size();
      nTracksAfter += nAfter;
      if (merge && nAfter < nin)
         nMergeEvents++;
      nFitted += te->GetFittedTracks().size();
   }
   std::cout << "[" << tag << "] events=" << n << "  multi-track events=" << nMultiTrackEvents
             << "  input tracks=" << nInTracks << "  tracks after=" << nTracksAfter << "  fitted=" << nFitted
             << (merge ? ("  merge-events=" + std::to_string(nMergeEvents)) : "") << "\n";
   f->Close();
}

void testMerge_a1975(TString runName = "run_0106", Long64_t nEv = 150, TString recoDir = "/mnt/f/a1975/reco/")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   TString dir = getenv("VMCWORKDIR");
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   TString geoManFile = dir + "/geometry/ATTPC_H1bar_geomanager.root";
   if (gROOT->FindObject("FAIRGeom") == nullptr) {
      TFile *gf = TFile::Open(geoManFile);
      gf->Get("FAIRGeom");
   }
   TString recoFile = recoDir + runName + "_reco.root";

   int a, b, c, d, e;
   runPass("merge OFF  ", 0, recoFile, nEv, a, b, c, d, e);
   int a2, b2, c2, d2, e2;
   runPass("merge ON   ", 1, recoFile, nEv, a2, b2, c2, d2, e2);
   int a3, b3, c3, d3, e3;
   runPass("force-merge", 2, recoFile, nEv, a3, b3, c3, d3, e3);
   std::cout << "\nDefault criteria  : tracks " << a << " -> " << d2 << ", fitted " << b << " -> " << b2
             << ", merge-events " << c2 << "\n";
   std::cout << "Force (loose) test: " << e << " multi-track events -> collapsed in " << c3
             << " events (tracks " << a << " -> " << d3 << "); proves merge mechanics fire.\n";
}
