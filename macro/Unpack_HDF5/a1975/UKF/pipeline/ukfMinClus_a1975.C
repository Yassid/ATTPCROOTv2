/// @file ukfMinClus_a1975.C
/// @brief Measure UKF fitted-track yield vs MinClusters on a1975 (p,p). The viewer
/// showed MinClusters(10) drops short tracks that genfit (MinClusters 4) keeps. Scan
/// the threshold and report yield + a crude quality proxy (fraction with physical KE
/// 1-60 MeV and median KE) so we can recover yield without admitting junk.
///   root -l -b -q 'ukfMinClus_a1975.C("run_0106", 3000, "/mnt/f/a1975/reco/")'

double med(std::vector<double> v) { if (v.empty()) return -1; std::sort(v.begin(), v.end()); return v[v.size()/2]; }

void onePass(int minClus, TTree *t, TClonesArray *peArr, Long64_t n)
{
   const double kE_C = 1.602176634e-19;
   auto eloss = std::make_unique<AtTools::AtELossCATIMA>(9.0e-5);
   eloss->SetProjectile(1, 1, 938.27208816 / 931.49410242);
   std::vector<std::tuple<int,int,int>> mat; mat.push_back({1,1,1}); eloss->SetMaterial(mat);
   auto fit = new EventFit::AtFitterUKF(kE_C, 938.27208816, std::move(eloss));
   fit->SetBField(ROOT::Math::XYZVector(0, 0, -2.85));
   fit->SetUKFParameters(1e-3, 2.0, 0.0);
   fit->SetMeasurementSigma(2.0);
   fit->SetMomentumSigmaFrac(0.3);
   fit->SetEnableEnergyStraggling(false);
   fit->SetMinClusters(minClus);
   fit->SetZPadPlane(1000.0);

   std::vector<double> kes; int nFit = 0, nPhys = 0;
   for (Long64_t i = 0; i < n; ++i) {
      t->GetEntry(i);
      auto *pe = (AtPatternEvent *)peArr->At(0);
      if (!pe) continue;
      auto te = std::make_unique<AtTrackingEvent>();
      fit->FitEvent(te.get(), pe);
      for (const auto &ft : te->GetFittedTracks()) {
         nFit++;
         double ke = ft->GetKinematics().kineticEnergy;
         kes.push_back(ke);
         if (ke > 1.0 && ke < 60.0) nPhys++;
      }
   }
   printf("  MinClusters=%-2d : fitted=%4d  physical-KE(1-60)=%4d (%.0f%%)  median KE=%.2f MeV\n", minClus, nFit, nPhys,
          nFit ? 100.0 * nPhys / nFit : 0, med(kes));
   delete fit;
}

void ukfMinClus_a1975(TString runName = "run_0106", Long64_t nEv = 3000, TString recoDir = "/mnt/f/a1975/reco/")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   auto *f = TFile::Open(recoDir + runName + "_reco.root");
   auto *t = (TTree *)f->Get("cbmsim");
   auto *peArr = new TClonesArray("AtPatternEvent");
   t->SetBranchAddress("AtPatternEvent", &peArr);
   Long64_t n = std::min(nEv, t->GetEntries());
   printf("\n=== UKF yield vs MinClusters, %s, %lld events (genfit uses 4) ===\n", runName.Data(), n);
   for (int mc : {10, 7, 5, 4, 3}) onePass(mc, t, peArr, n);
}
