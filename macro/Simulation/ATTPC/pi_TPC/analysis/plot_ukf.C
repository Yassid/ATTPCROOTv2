/// @file plot_ukf.C
/// @brief Build a 6-panel diagnostic for the pi-in-TPC UKF run.
///
/// Reads ./data/attpcsim.root (MC truth) and ./data/output_ukf_only.root
/// (fitted tracks) and produces ./data/ukf_diagnostics.png with:
///   1. KE residual (fit - truth, MeV) + Gaussian fit
///   2. KE_fit vs KE_truth scatter + diagonal
///   3. theta residual (deg) + Gaussian fit
///   4. fitted vertex (x, y) with truth marker
///   5. fitted vertex z
///   6. reduced chi^2 distribution
///
/// Run: root -b -q plot_ukf.C

void plot_ukf()
{
   gStyle->SetOptStat(0);
   gStyle->SetOptFit(1111);
   gStyle->SetTitleSize(0.05, "XYZ");
   gStyle->SetLabelSize(0.045, "XYZ");
   gStyle->SetTitleSize(0.055, "T");

   TFile fSim("data/attpcsim.root");
   TFile fUKF("data/output_ukf_only.root");
   auto *tSim = (TTree *)fSim.Get("cbmsim");
   auto *tUKF = (TTree *)fUKF.Get("cbmsim");
   if (!tSim || !tUKF) {
      std::cout << "missing trees\n";
      return;
   }

   TClonesArray *trks = new TClonesArray("AtMCTrack");
   TClonesArray *teArr = new TClonesArray("AtTrackingEvent");
   tSim->SetBranchAddress("MCTrack", &trks);
   tUKF->SetBranchAddress("AtTrackingEvent", &teArr);

   // Histos
   TH1F *hKEres = new TH1F("hKEres", "KE residual;KE_{fit} - KE_{MC} (MeV);counts", 80, -25., 25.);
   TH2F *hKEcorr = new TH2F("hKEcorr", "KE_{fit} vs KE_{MC};KE_{MC} (MeV);KE_{fit} (MeV)", 60, 0., 55., 60, 0., 55.);
   TH1F *hThres = new TH1F("hThres", "#theta residual;#theta_{fit} - #theta_{MC} (deg);counts", 80, -8., 8.);
   TH2F *hVxy = new TH2F("hVxy", "Reconstructed vertex (x,y);x (mm);y (mm)", 80, -100., 100., 80, -100., 100.);
   TH1F *hVz = new TH1F("hVz", "Reconstructed vertex z;z (mm);counts", 80, 0., 1000.);
   TH1F *hChi2 = new TH1F("hChi2", "Reduced #chi^{2};#chi^{2}/ndf;counts", 100, 0., 0.5);

   Long64_t n = std::min(tSim->GetEntries(), tUKF->GetEntries());
   for (Long64_t i = 0; i < n; ++i) {
      tSim->GetEntry(i);
      tUKF->GetEntry(i);
      if (trks->GetEntries() == 0) continue;
      auto *mc = (AtMCTrack *)trks->At(0);
      if (mc->GetPdgCode() != 211 && mc->GetPdgCode() != -211) continue;
      double KEmc = (mc->GetEnergy() - mc->GetMass()) * 1000.;
      double pmc = std::sqrt(mc->GetPx() * mc->GetPx() + mc->GetPy() * mc->GetPy() + mc->GetPz() * mc->GetPz());
      double thMC = std::acos(mc->GetPz() / pmc) * 180. / M_PI;

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
      if (!best || !best->GetTrackMetadata()) continue;

      double ndf = best->GetTrackMetadata()->GetNdf();
      double redChi = ndf > 0 ? best->GetTrackMetadata()->GetChi2() / ndf : 0;

      const auto &k = best->GetKinematics();
      double KEfit = k.kineticEnergy; // MeV
      double thFit = k.theta * 180. / M_PI;
      const auto &v = best->GetVertex();

      hKEres->Fill(KEfit - KEmc);
      hKEcorr->Fill(KEmc, KEfit);
      hThres->Fill(thFit - thMC);
      hVxy->Fill(v.X(), v.Y());
      hVz->Fill(v.Z());
      hChi2->Fill(redChi);
   }

   auto *c = new TCanvas("c", "UKF pi-in-TPC diagnostics", 1500, 1000);
   c->Divide(3, 2);

   c->cd(1);
   gPad->SetGrid();
   hKEres->SetFillColorAlpha(kAzure - 4, 0.5);
   hKEres->SetLineColor(kAzure + 2);
   hKEres->Fit("gaus", "Q");
   if (auto *f = hKEres->GetFunction("gaus")) f->SetLineColor(kRed + 1);
   hKEres->Draw("hist");
   if (auto *f = hKEres->GetFunction("gaus")) f->Draw("same");

   c->cd(2);
   gPad->SetGrid();
   hKEcorr->SetMarkerStyle(20);
   hKEcorr->SetMarkerSize(0.5);
   hKEcorr->SetMarkerColor(kAzure + 2);
   hKEcorr->Draw("colz");
   TLine *diag = new TLine(0, 0, 55, 55);
   diag->SetLineColor(kRed + 1);
   diag->SetLineWidth(2);
   diag->SetLineStyle(2);
   diag->Draw("same");

   c->cd(3);
   gPad->SetGrid();
   hThres->SetFillColorAlpha(kSpring - 7, 0.5);
   hThres->SetLineColor(kGreen + 2);
   hThres->Fit("gaus", "Q");
   if (auto *f = hThres->GetFunction("gaus")) f->SetLineColor(kRed + 1);
   hThres->Draw("hist");
   if (auto *f = hThres->GetFunction("gaus")) f->Draw("same");

   c->cd(4);
   gPad->SetGrid();
   hVxy->Draw("colz");
   TMarker *truthVxy = new TMarker(0., 0., 29);
   truthVxy->SetMarkerColor(kRed + 1);
   truthVxy->SetMarkerSize(2.0);
   truthVxy->Draw("same");

   c->cd(5);
   gPad->SetGrid();
   hVz->SetFillColorAlpha(kOrange + 1, 0.5);
   hVz->SetLineColor(kOrange + 3);
   hVz->Draw("hist");
   TLine *truthVz = new TLine(500., 0., 500., hVz->GetMaximum());
   truthVz->SetLineColor(kRed + 1);
   truthVz->SetLineWidth(2);
   truthVz->SetLineStyle(2);
   truthVz->Draw("same");

   c->cd(6);
   gPad->SetGrid();
   hChi2->SetFillColorAlpha(kViolet - 7, 0.5);
   hChi2->SetLineColor(kViolet + 2);
   hChi2->Draw("hist");

   c->SaveAs("data/ukf_diagnostics.png");
   c->SaveAs("data/ukf_diagnostics.pdf");
   std::cout << "Wrote data/ukf_diagnostics.png and .pdf\n";
}
