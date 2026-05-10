/// @brief 16C(p,p) recoil-proton kinematics overlay — MC truth vs UKF
/// reconstructed in the (θ_lab, KE) plane. The reaction is two-body
/// elastic; truth points should fall on a thin curve and the fit should
/// reproduce it.
///
/// Run from 16C_pp/: `root -b -q make_kinematics.C`

void make_kinematics()
{
   gStyle->SetOptStat(0);
   gStyle->SetTitleSize(0.055, "T");
   gStyle->SetTitleSize(0.05, "XY");
   gStyle->SetLabelSize(0.045, "XY");
   gStyle->SetPadLeftMargin(0.13);
   gStyle->SetPadBottomMargin(0.13);

   TFile fSim("data/attpcsim.root");
   TFile fUKF("data/output_ukf_perf.root");
   auto *tSim = (TTree *)fSim.Get("cbmsim");
   auto *tUKF = (TTree *)fUKF.Get("cbmsim");
   if (!tSim || !tUKF) { std::cerr << "missing trees\n"; return; }

   TClonesArray *trks = new TClonesArray("AtMCTrack");
   TClonesArray *teArr = new TClonesArray("AtTrackingEvent");
   tSim->SetBranchAddress("MCTrack", &trks);
   tUKF->SetBranchAddress("AtTrackingEvent", &teArr);

   TH2F *hTruth = new TH2F("hTruth", "MC truth;#theta_{lab} (deg);KE (MeV)", 90, 30., 90., 100, 0., 25.);
   TH2F *hReco = new TH2F("hReco", "UKF reconstructed;#theta_{lab} (deg);KE (MeV)", 90, 30., 90., 100, 0., 25.);
   TH2F *hOver = new TH2F("hOver", "Truth (black) vs UKF (red);#theta_{lab} (deg);KE (MeV)",
                          90, 30., 90., 100, 0., 25.);

   Long64_t n = std::min(tSim->GetEntries(), tUKF->GetEntries());
   int nProton = 0, nFit = 0;
   for (Long64_t i = 0; i < n; ++i) {
      tSim->GetEntry(i);
      tUKF->GetEntry(i);
      AtMCTrack *p = nullptr;
      for (int j = 0; j < trks->GetEntries(); ++j) {
         auto *mc = (AtMCTrack *)trks->At(j);
         if (mc->GetPdgCode() == 2212 && (mc->GetEnergy() - mc->GetMass()) > 0.001) { p = mc; break; }
      }
      if (!p) continue;
      double pmc_GeV = std::sqrt(p->GetPx() * p->GetPx() + p->GetPy() * p->GetPy() + p->GetPz() * p->GetPz());
      double KEmc = (p->GetEnergy() - p->GetMass()) * 1000.;
      double thMC = std::acos(p->GetPz() / pmc_GeV) * 180. / M_PI;
      hTruth->Fill(thMC, KEmc);
      ++nProton;

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
         if (chi < bestChi) { bestChi = chi; best = t.get(); }
      }
      if (!best) continue;
      if (bestChi > 10.0) continue;
      double KEfit = best->GetKinematics().kineticEnergy;
      double thFit = best->GetKinematics().theta * 180. / M_PI;
      if (KEfit < 0.5 || KEfit > 50.0) continue;
      hReco->Fill(thFit, KEfit);
      ++nFit;
   }

   auto *c = new TCanvas("c", "16C(p,p) kinematics", 1500, 600);
   c->Divide(2, 1, 0.005, 0.02);

   c->cd(1);
   gPad->SetGrid();
   hTruth->SetMarkerStyle(20);
   hTruth->SetMarkerSize(0.4);
   hTruth->SetMarkerColor(kBlack);
   hTruth->Draw();

   c->cd(2);
   gPad->SetGrid();
   hReco->SetMarkerStyle(20);
   hReco->SetMarkerSize(0.4);
   hReco->SetMarkerColor(kRed + 1);
   hReco->Draw();

   c->SaveAs("data/kinematics_perf.png");
   c->SaveAs("data/kinematics_perf.pdf");

   // Overlay
   auto *cOv = new TCanvas("cOv", "16C(p,p) kinematics overlay", 900, 700);
   gPad->SetGrid();
   hTruth->SetMarkerSize(0.5);
   hTruth->SetMarkerColor(kBlack);
   hTruth->SetTitle("Truth (black) vs UKF reconstructed (red);#theta_{lab} (deg);KE (MeV)");
   hTruth->Draw();
   hReco->SetMarkerSize(0.4);
   hReco->SetMarkerColor(kRed + 1);
   hReco->Draw("same");
   auto *leg = new TLegend(0.55, 0.74, 0.88, 0.88);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->AddEntry(hTruth, Form("MC truth proton (%d)", nProton), "p");
   leg->AddEntry(hReco, Form("UKF reconstructed (%d)", nFit), "p");
   leg->Draw();
   cOv->SaveAs("data/kinematics_overlay_perf.png");
   cOv->SaveAs("data/kinematics_overlay_perf.pdf");

   std::cout << "Wrote data/kinematics_perf.{png,pdf} and kinematics_overlay_perf.{png,pdf}\n";
   std::cout << "MC protons: " << nProton << "  Reconstructed: " << nFit
             << "  (" << 100. * nFit / std::max(1, nProton) << "%)\n";
}
