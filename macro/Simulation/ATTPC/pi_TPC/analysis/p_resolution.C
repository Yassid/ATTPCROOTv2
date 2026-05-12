/// @file p_resolution.C
/// @brief Print σ_p/p per θ bin for the high-energy regime where σ_KE/KE
/// is not a useful proxy. Computes p_reco from KE_reco + mass and compares
/// to MC truth p.
///
/// Usage:
///   root -b -q 'analysis/p_resolution.C("data/output_ukf_only_hipE_B2T.root",
///       "data/attpcsim_hipE_B2T.root","B=2 T, 800 MeV/c")'

void p_resolution(const char *ukfFile = "data/output_ukf_only.root",
                  const char *simFile = "data/attpcsim.root",
                  const char *label = "")
{
   TFile fSim(simFile);
   TFile fUKF(ukfFile);
   auto *tSim = (TTree *)fSim.Get("cbmsim");
   auto *tUKF = (TTree *)fUKF.Get("cbmsim");

   auto *trks = new TClonesArray("AtMCTrack");
   auto *te = new TClonesArray("AtTrackingEvent");
   tSim->SetBranchAddress("MCTrack", &trks);
   tUKF->SetBranchAddress("AtTrackingEvent", &te);

   const int NB = 9;
   double edges[NB + 1] = {5., 25., 45., 65., 80., 100., 115., 135., 155., 175.};
   std::vector<TH1F *> hDP(NB), hP(NB);
   std::vector<int> nThr(NB, 0);
   for (int b = 0; b < NB; ++b) {
      hDP[b] = new TH1F(Form("hDP_%d", b), "", 100, -1.5, 1.5); // Δp/p
      hP[b] = new TH1F(Form("hP_%d", b), "", 100, 0., 2000.);
   }

   const double mass_pi = 139.57039; // MeV/c^2
   Long64_t n = std::min(tSim->GetEntries(), tUKF->GetEntries());
   for (Long64_t i = 0; i < n; ++i) {
      tSim->GetEntry(i);
      tUKF->GetEntry(i);
      if (trks->GetEntries() == 0) continue;
      auto *mc = (AtMCTrack *)trks->At(0);
      if (std::abs(mc->GetPdgCode()) != 211) continue;
      double pmc = std::sqrt(mc->GetPx() * mc->GetPx() + mc->GetPy() * mc->GetPy() + mc->GetPz() * mc->GetPz()) * 1000.;
      double thMC = std::acos(mc->GetPz() / (pmc / 1000.)) * 180. / M_PI;
      int b = -1;
      for (int k = 0; k < NB; ++k)
         if (thMC >= edges[k] && thMC < edges[k + 1]) { b = k; break; }
      if (b < 0) continue;
      ++nThr[b];
      hP[b]->Fill(pmc);
      if (te->GetEntries() == 0) continue;
      auto *trkEvt = (AtTrackingEvent *)te->At(0);
      auto &fitted = trkEvt->GetFittedTracks();
      AtFittedTrack *best = nullptr;
      double bestChi = 1e30;
      for (auto &t : fitted) {
         if (!t->GetTrackMetadata()) continue;
         double ndf = t->GetTrackMetadata()->GetNdf();
         double c = ndf > 0 ? t->GetTrackMetadata()->GetChi2() / ndf : 1e30;
         if (c < bestChi) { bestChi = c; best = t.get(); }
      }
      if (!best) continue;
      double KEfit = best->GetKinematics().kineticEnergy;
      double Efit = KEfit + mass_pi;
      double pfit = std::sqrt(Efit * Efit - mass_pi * mass_pi);
      hDP[b]->Fill((pfit - pmc) / pmc);
   }

   std::cout << "\n=== p_resolution [" << label << "] file=" << ukfFile << " ===\n";
   std::cout << "theta(deg)   Nthr  Nfit   <p_mc>      bias       σ_p/p\n";
   std::cout << std::string(70, '-') << "\n";
   for (int b = 0; b < NB; ++b) {
      double meanP = hP[b]->GetEntries() > 0 ? hP[b]->GetMean() : 0;
      double mean = hDP[b]->GetMean(), rms = hDP[b]->GetRMS();
      if (hDP[b]->GetEntries() >= 15) {
         hDP[b]->Fit("gaus", "Q0", "", mean - 2.0 * rms, mean + 2.0 * rms);
         auto *f = hDP[b]->GetFunction("gaus");
         if (f) { mean = f->GetParameter(1); rms = f->GetParameter(2); }
      }
      printf("%4.0f-%-4.0f    %4d  %4.0f   %7.1f    %+6.2f     %5.2f (%5.1f%%)\n",
             edges[b], edges[b + 1], nThr[b], hDP[b]->GetEntries(), meanP, mean, rms, rms * 100.);
   }
}
