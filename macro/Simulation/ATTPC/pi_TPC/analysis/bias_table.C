/// @file bias_table.C
/// @brief Print KE bias and σ/E per theta bin for a given UKF output file.
///
/// Usage:
///   root -b -q 'analysis/bias_table.C("data/output_ukf_only_mom03.root","mom=0.3")'

void bias_table(const char *ukfFile = "data/output_ukf_only.root",
                const char *label = "baseline",
                const char *simFile = "data/attpcsim.root")
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
   std::vector<TH1F *> hKE(NB), hKEcl(NB);
   std::vector<int> nThrown(NB, 0);
   std::vector<double> sumKE(NB, 0.), sumKE2(NB, 0.);
   for (int b = 0; b < NB; ++b) {
      hKE[b] = new TH1F(Form("hKE_%d", b), "", 100, -30., 30.);
      hKEcl[b] = new TH1F(Form("hKEcl_%d", b), "", 100, -30., 30.);
   }

   Long64_t n = std::min(tSim->GetEntries(), tUKF->GetEntries());
   for (Long64_t i = 0; i < n; ++i) {
      tSim->GetEntry(i);
      tUKF->GetEntry(i);
      if (trks->GetEntries() == 0) continue;
      auto *mc = (AtMCTrack *)trks->At(0);
      if (std::abs(mc->GetPdgCode()) != 211) continue;
      double KEmc = (mc->GetEnergy() - mc->GetMass()) * 1000.;
      double pmc = std::sqrt(mc->GetPx() * mc->GetPx() + mc->GetPy() * mc->GetPy() + mc->GetPz() * mc->GetPz());
      double thMC = std::acos(mc->GetPz() / pmc) * 180. / M_PI;
      int b = -1;
      for (int k = 0; k < NB; ++k)
         if (thMC >= edges[k] && thMC < edges[k + 1]) { b = k; break; }
      if (b < 0) continue;
      ++nThrown[b];
      sumKE[b] += KEmc;
      sumKE2[b] += KEmc * KEmc;
      if (te->GetEntries() == 0) continue;
      auto *trkEvt = (AtTrackingEvent *)te->At(0);
      auto &fitted = trkEvt->GetFittedTracks();
      if (fitted.empty()) continue;
      AtFittedTrack *best = nullptr;
      double bestChi = 1e30;
      for (auto &t : fitted) {
         if (!t->GetTrackMetadata()) continue;
         double ndf = t->GetTrackMetadata()->GetNdf();
         double chi = ndf > 0 ? t->GetTrackMetadata()->GetChi2() / ndf : 1e30;
         if (chi < bestChi) { bestChi = chi; best = t.get(); }
      }
      if (!best) continue;
      double KEfit = best->GetKinematics().kineticEnergy;
      double KEcl = best->GetKinematicsXtr().kineticEnergy;
      hKE[b]->Fill(KEfit - KEmc);
      hKEcl[b]->Fill(KEcl - KEmc);
   }

   std::cout << "\n=== bias_table [" << label << "] file=" << ukfFile << " ===\n";
   std::cout << "theta(deg)  Nthr  Nfit   <KEmc>   bias_vtx  sigma_vtx  bias_clu  sigma_clu\n";
   std::cout << std::string(86, '-') << "\n";
   for (int b = 0; b < NB; ++b) {
      if (hKE[b]->GetEntries() < 15) {
         printf("%4.0f-%-4.0f   %4d  %4.0f   ----     ----      ----       ----      ----\n",
                edges[b], edges[b + 1], nThrown[b], hKE[b]->GetEntries());
         continue;
      }
      double meanKE = sumKE[b] / nThrown[b];
      double mean = hKE[b]->GetMean();
      double rms = hKE[b]->GetRMS();
      hKE[b]->Fit("gaus", "Q0", "", mean - 2.5 * rms, mean + 2.5 * rms);
      auto *f = hKE[b]->GetFunction("gaus");
      double bias = f ? f->GetParameter(1) : mean;
      double sigma = f ? f->GetParameter(2) : rms;
      double meanc = hKEcl[b]->GetMean();
      double rmsc = hKEcl[b]->GetRMS();
      hKEcl[b]->Fit("gaus", "Q0", "", meanc - 2.5 * rmsc, meanc + 2.5 * rmsc);
      auto *fc = hKEcl[b]->GetFunction("gaus");
      double biasc = fc ? fc->GetParameter(1) : meanc;
      double sigmac = fc ? fc->GetParameter(2) : rmsc;
      printf("%4.0f-%-4.0f   %4d  %4.0f   %5.1f    %+5.2f     %5.2f      %+5.2f     %5.2f\n",
             edges[b], edges[b + 1], nThrown[b], hKE[b]->GetEntries(), meanKE,
             bias, sigma, biasc, sigmac);
   }
}
