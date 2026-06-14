/// @file mergeQuality_a1975.C
/// @brief Quantify the EFFECT of continuity merging: for events where merging changes
/// the track set, compare fragment fits (merge OFF) vs the merged fit (merge ON), and
/// classify each merge by the QUALITY of the result (not a naive track count):
///   HARMFUL  : a track that fit well (goodFit) as fragments fits BADLY when merged
///              (chi2/ndf > cut or no convergence) — a wrong join.
///   DEDUP    : >=2 good fragment-fits collapse to ONE still-good merged fit — the
///              fragments were one PRA-split track, correctly rejoined.
///   RECOVER  : fragments produced NO good fit, the merged track does — a recovery.
///   NEUTRAL  : everything else.
/// proton channel:  root -l -b -q 'mergeQuality_a1975.C("run_0106",5000)'
/// deuteron (p,d):  root -l -b -q 'mergeQuality_a1975.C("run_0106",20000,"/mnt/f/a1975/reco/",1000010020,2.0135532,1,"deuteron_H2_catima.txt","pid/deuteron_band.json")'

struct TrkQ {
   int nGood = 0;
   double bestChi2ndf = 1e9; // smallest chi2/ndf among converged good fits
};

TrkQ measure(EventFit::AtGenfitter *fitter, AtPatternEvent *pe, AtTrackingEvent *teOut = nullptr)
{
   TrkQ q;
   auto te = std::make_unique<AtTrackingEvent>();
   AtTrackingEvent *use = teOut ? teOut : te.get();
   if (!teOut) fitter->FitEvent(use, pe);
   for (const auto &ft : use->GetFittedTracks()) {
      if (const auto &m = ft->GetTrackMetadata()) {
         double ndf = m->GetNdf();
         double c = ndf > 0 ? m->GetChi2() / ndf : 1e9;
         if (m->GetGoodFit()) { q.nGood++; q.bestChi2ndf = std::min(q.bestChi2ndf, c); }
      }
   }
   return q;
}

void mergeQuality_a1975(TString runName = "run_0106", Long64_t nEv = 5000, TString recoDir = "/mnt/f/a1975/reco/",
                        Int_t pdg = 2212, Double_t massAmu = 938.27208816 / 931.49410242, Int_t Z = 1,
                        TString elossName = "", TString pidGate = "")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   TString dir = getenv("VMCWORKDIR");
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   if (gROOT->FindObject("FAIRGeom") == nullptr)
      TFile::Open(dir + "/geometry/ATTPC_H1bar_geomanager.root")->Get("FAIRGeom");
   std::string elossFile = elossName.Length() ? ((std::string)dir.Data() + "/resources/energy_loss/" + elossName.Data()) : std::string("");

   auto *f = TFile::Open(recoDir + runName + "_reco.root");
   auto *t = (TTree *)f->Get("cbmsim");
   auto *peArr = new TClonesArray("AtPatternEvent");
   t->SetBranchAddress("AtPatternEvent", &peArr);

   auto mk = [&](bool merge) {
      auto *fit = new EventFit::AtGenfitter(-2.85, pdg, massAmu, Z, elossFile, kTRUE, 2, 5);
      fit->SetZPadPlane(1000.0);
      fit->SetMeasSigma(4.0);
      if (pidGate.Length() && !gSystem->AccessPathName(pidGate.Data())) fit->SetPIDGate(pidGate.Data());
      if (merge) fit->SetMergeContinuity(kTRUE); // tuned defaults (centreDist 15mm)
      fit->Init();
      return fit;
   };
   auto *fOff = mk(false), *fOn = mk(true);

   Long64_t n = std::min(nEv, t->GetEntries());
   int mergeEvents = 0, harmful = 0, dedup = 0, recover = 0, neutral = 0;
   for (Long64_t i = 0; i < n; ++i) {
      t->GetEntry(i);
      auto *pe = (AtPatternEvent *)peArr->At(0);
      if (!pe || pe->GetTrackCand().size() < 2) continue;

      auto teOn = std::make_unique<AtTrackingEvent>();
      fOn->FitEvent(teOn.get(), pe);
      if ((int)teOn->GetTrackArray().size() >= (int)pe->GetTrackCand().size()) continue; // no merge fired
      mergeEvents++;

      TrkQ off = measure(fOff, pe);
      TrkQ on = measure(nullptr, pe, teOn.get());

      const double cut = 5.0;
      TString cls;
      if (off.nGood == 0 && on.nGood > 0) { recover++; cls = "RECOVER"; }
      else if (off.nGood > 0 && (on.nGood == 0 || on.bestChi2ndf > cut)) { harmful++; cls = "HARMFUL"; }
      else if (off.nGood >= 2 && on.nGood >= 1 && on.bestChi2ndf <= cut) { dedup++; cls = "DEDUP"; }
      else { neutral++; cls = "NEUTRAL"; }

      std::cout << "  ev " << i << " [" << cls << "] OFF good=" << off.nGood << " bestChi2ndf=" << off.bestChi2ndf
                << "  ->  ON good=" << on.nGood << " bestChi2ndf=" << on.bestChi2ndf << "\n";
   }

   std::cout << "\n=== " << runName << ", " << n << " events, pdg=" << pdg
             << (pidGate.Length() ? (" gate=" + pidGate) : " (ungated)") << " ===\n";
   std::cout << "merge fired   : " << mergeEvents << " events\n";
   std::cout << "  RECOVER     : " << recover << "  (no good fragment fit -> good merged fit)\n";
   std::cout << "  DEDUP       : " << dedup << "  (split track correctly rejoined, still good)\n";
   std::cout << "  HARMFUL     : " << harmful << "  (good fragment -> bad merged fit)\n";
   std::cout << "  NEUTRAL     : " << neutral << "\n";
   std::cout << "NET fit quality: " << (harmful == 0 ? "no degradation" : (std::to_string(harmful) + " degraded"))
             << "; " << (recover + dedup) << " events improved/deduplicated\n";
}
