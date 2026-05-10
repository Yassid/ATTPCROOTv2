/// @brief Quantify KE resolution overall and in KE bins.
void ke_resolution()
{
   TFile fSim("data/attpcsim.root");
   TFile fUKF("data/output_ukf_only.root");
   auto *tSim = (TTree *)fSim.Get("cbmsim");
   auto *tUKF = (TTree *)fUKF.Get("cbmsim");

   TClonesArray *trks = new TClonesArray("AtMCTrack");
   TClonesArray *teArr = new TClonesArray("AtTrackingEvent");
   tSim->SetBranchAddress("MCTrack", &trks);
   tUKF->SetBranchAddress("AtTrackingEvent", &teArr);

   const int NB = 5;
   double edges[NB + 1] = {5., 14., 23., 32., 41., 50.};
   std::vector<TH1F *> hRes(NB);
   for (int b = 0; b < NB; ++b) {
      hRes[b] = new TH1F(Form("hRes_%d", b),
                         Form("KE %.0f-%.0f;dKE (MeV);counts", edges[b], edges[b + 1]), 80, -25., 25.);
   }
   TH1F *hAll = new TH1F("hAll", "all KE;dKE (MeV);counts", 100, -25., 25.);

   Long64_t n = std::min(tSim->GetEntries(), tUKF->GetEntries());
   int nFit = 0;
   for (Long64_t i = 0; i < n; ++i) {
      tSim->GetEntry(i);
      tUKF->GetEntry(i);
      if (trks->GetEntries() == 0) continue;
      auto *mc = (AtMCTrack *)trks->At(0);
      if (std::abs(mc->GetPdgCode()) != 211) continue;
      double KEmc = (mc->GetEnergy() - mc->GetMass()) * 1000.;

      if (teArr->GetEntries() == 0) continue;
      auto *te = (AtTrackingEvent *)teArr->At(0);
      auto &fitted = te->GetFittedTracks();
      if (fitted.empty()) continue;

      AtFittedTrack *best = nullptr;
      double bestChi = 1e30;
      for (auto &t : fitted) {
         if (!t->GetTrackMetadata()) continue;
         double ndf = t->GetTrackMetadata()->GetNdf();
         double chi = ndf > 0 ? t->GetTrackMetadata()->GetChi2() / ndf : 1e30;
         if (chi < bestChi) {
            bestChi = chi;
            best = t.get();
         }
      }
      if (!best) continue;
      double KEfit = best->GetKinematics().kineticEnergy;
      double dKE = KEfit - KEmc;
      hAll->Fill(dKE);
      for (int b = 0; b < NB; ++b) {
         if (KEmc >= edges[b] && KEmc < edges[b + 1]) {
            hRes[b]->Fill(dKE);
            break;
         }
      }
      ++nFit;
   }

   std::cout << "\n=== KE resolution (theta=45 fixed, B=0.5T, P10@1bar, KE 5-50 MeV pi+) ===\n";
   std::cout << "Fitted entries: " << nFit << " / " << n << "\n\n";

   auto fitGaus = [](TH1F *h, const char *label, double Emid) {
      if (h->GetEntries() < 20) {
         std::cout << label << " (entries=" << h->GetEntries() << ", skip)\n";
         return;
      }
      double rms = h->GetRMS();
      double mean = h->GetMean();
      double lo = std::max(-20., mean - 2.5 * rms);
      double hi = std::min(20., mean + 2.5 * rms);
      h->Fit("gaus", "Q0", "", lo, hi);
      auto *f = h->GetFunction("gaus");
      double mu = f ? f->GetParameter(1) : mean;
      double sg = f ? f->GetParameter(2) : rms;
      double rel = Emid > 0 ? sg / Emid * 100. : 0;
      std::cout << label << "  N=" << std::setw(4) << (int)h->GetEntries()
                << "  mean=" << std::fixed << std::setprecision(2) << std::setw(6) << mean
                << "  RMS=" << std::setw(5) << rms
                << "  Gaus mu=" << std::setw(6) << mu << " sigma=" << std::setw(5) << sg;
      if (Emid > 0) std::cout << "   sigma/E = " << std::setw(5) << std::setprecision(1) << rel << "%";
      std::cout << "\n";
   };

   fitGaus(hAll, "All KE          ", 27.5);
   for (int b = 0; b < NB; ++b) {
      double mid = 0.5 * (edges[b] + edges[b + 1]);
      char buf[64];
      snprintf(buf, sizeof(buf), "KE %4.1f-%4.1f MeV", edges[b], edges[b + 1]);
      fitGaus(hRes[b], buf, mid);
   }
}
