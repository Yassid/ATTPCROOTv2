/// @file ukfMSscan_a1975.C
/// @brief Measure the effect of UKF multiple scattering on a1975 (p,p): fit the same
/// proton events with MS off vs on (Highland, H2 X0) and compare yield + median KE.
/// Prediction: negligible in H2 1 bar (per-step MS variance ~1e-9 rad^2 << the
/// fAngModelNoise floor 1e-6). This confirms it quantitatively.
///   root -l -b -q 'ukfMSscan_a1975.C("run_0106", 2000, "/mnt/f/a1975/reco/")'

double med(std::vector<double> v) { if (v.empty()) return -1; std::sort(v.begin(), v.end()); return v[v.size()/2]; }

void onePass(bool ms, TTree *t, TClonesArray *peArr, Long64_t n)
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
   fit->SetMinClusters(10);
   fit->SetZPadPlane(1000.0);
   if (ms) { fit->SetEnableMultipleScattering(true); fit->SetRadiationLength(7.6e6); }

   std::vector<double> kes; int nFit = 0;
   for (Long64_t i = 0; i < n; ++i) {
      t->GetEntry(i);
      auto *pe = (AtPatternEvent *)peArr->At(0);
      if (!pe) continue;
      auto te = std::make_unique<AtTrackingEvent>();
      fit->FitEvent(te.get(), pe);
      for (const auto &ft : te->GetFittedTracks()) { nFit++; kes.push_back(ft->GetKinematics().kineticEnergy); }
   }
   printf("  UKF MS %-3s : fitted=%4d  median KE=%.3f MeV\n", ms ? "ON" : "OFF", nFit, med(kes));
   delete fit;
}

void ukfMSscan_a1975(TString runName = "run_0106", Long64_t nEv = 2000, TString recoDir = "/mnt/f/a1975/reco/")
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   auto *f = TFile::Open(recoDir + runName + "_reco.root");
   auto *t = (TTree *)f->Get("cbmsim");
   auto *peArr = new TClonesArray("AtPatternEvent");
   t->SetBranchAddress("AtPatternEvent", &peArr);
   Long64_t n = std::min(nEv, t->GetEntries());
   printf("\n=== UKF multiple-scattering effect, %s, %lld events ===\n", runName.Data(), n);
   onePass(false, t, peArr, n);
   onePass(true, t, peArr, n);
}
