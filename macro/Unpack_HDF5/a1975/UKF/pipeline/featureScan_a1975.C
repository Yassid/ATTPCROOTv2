/// @file featureScan_a1975.C
/// @brief Measure the EFFECT of the new genfit covariance knobs on a1975 (p,p):
/// for each (base sigma, diffusion) config, fit N proton-gated events and report the
/// median chi2/ndf and the goodFit rate. A well-calibrated covariance gives median
/// chi2/ndf ~ 1; the legacy flat 4 mm gives << 1 (errors over-estimated).
///   root -l -b -q 'featureScan_a1975.C("run_0106", 3000, "/mnt/f/a1975/reco/")'

double medianOf(std::vector<double> v)
{
   if (v.empty()) return -1;
   std::sort(v.begin(), v.end());
   return v[v.size() / 2];
}

void onePass(const char *tag, double sigma, double diffT, double diffL, TTree *t, TClonesArray *peArr, Long64_t n,
             TString pidGate)
{
   auto *fit = new EventFit::AtGenfitter(-2.85, 2212, 938.27208816 / 931.49410242, 1, "", kTRUE, 2, 5);
   fit->SetZPadPlane(1000.0);
   fit->SetMeasSigma(sigma);
   if (diffT > 0 || diffL > 0) fit->SetDiffusion(diffT, diffL);
   if (pidGate.Length() && !gSystem->AccessPathName(pidGate.Data())) fit->SetPIDGate(pidGate.Data());
   fit->Init();

   std::vector<double> chi2ndf;
   int nFit = 0, nGood = 0;
   for (Long64_t i = 0; i < n; ++i) {
      t->GetEntry(i);
      auto *pe = (AtPatternEvent *)peArr->At(0);
      if (!pe) continue;
      auto te = std::make_unique<AtTrackingEvent>();
      fit->FitEvent(te.get(), pe);
      for (const auto &ft : te->GetFittedTracks()) {
         nFit++;
         if (const auto &m = ft->GetTrackMetadata()) {
            double ndf = m->GetNdf();
            if (ndf > 0) chi2ndf.push_back(m->GetChi2() / ndf);
            if (m->GetGoodFit()) nGood++;
         }
      }
   }
   printf("  %-26s sigma=%.1f diffT=%.2f : fitted=%4d good=%4d (%.0f%%)  median chi2/ndf=%.3f\n", tag, sigma, diffT,
          nFit, nGood, nFit ? 100.0 * nGood / nFit : 0, medianOf(chi2ndf));
   delete fit;
}

void featureScan_a1975(TString runName = "run_0106", Long64_t nEv = 3000, TString recoDir = "/mnt/f/a1975/reco/",
                       TString pidGate = "pid/proton_band.json")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   TString dir = getenv("VMCWORKDIR");
   gSystem->Setenv("GEOMPATH", (dir + "/geometry/").Data());
   if (gROOT->FindObject("FAIRGeom") == nullptr)
      TFile::Open(dir + "/geometry/ATTPC_H1bar_geomanager.root")->Get("FAIRGeom");

   auto *f = TFile::Open(recoDir + runName + "_reco.root");
   auto *t = (TTree *)f->Get("cbmsim");
   auto *peArr = new TClonesArray("AtPatternEvent");
   t->SetBranchAddress("AtPatternEvent", &peArr);
   Long64_t n = std::min(nEv, t->GetEntries());

   printf("\n=== genfit measurement-covariance scan, %s, %lld events, gate=%s ===\n", runName.Data(), n,
          pidGate.Data());
   printf("(well-calibrated => median chi2/ndf ~ 1)\n");
   // base-sigma sweep (flat covariance)
   for (double s : {4.0, 3.0, 2.0, 1.5, 1.0})
      onePass("flat", s, 0.0, 0.0, t, peArr, n, pidGate);
   // diffusion model on top of a tighter base (transverse + longitudinal mm/sqrt(cm))
   onePass("diffusion(0.5,0.7)", 2.0, 0.5, 0.7, t, peArr, n, pidGate);
   onePass("diffusion(1.0,1.4)", 1.5, 1.0, 1.4, t, peArr, n, pidGate);
}
