/// @file make_kinematics.C
/// @brief pi-TPC kinematics overlay — MC truth vs UKF reconstructed in
/// the (θ_lab, KE) plane. Unlike 16C(p,p) where the proton sits on a
/// thin elastic-kinematics curve, the single-particle pi+ generator
/// fills a uniform box in θ ∈ [5°,175°] × KE ∈ [5,50] MeV. The
/// reconstructed scatter should fill the same box (acceptance) with
/// the appropriate spread (resolution).
///
/// Run from pi_TPC/: `root -b -q analysis/make_kinematics.C`

void make_kinematics()
{
   gStyle->SetOptStat(0);
   gStyle->SetTitleSize(0.055, "T");
   gStyle->SetTitleSize(0.05, "XY");
   gStyle->SetLabelSize(0.045, "XY");
   gStyle->SetPadLeftMargin(0.13);
   gStyle->SetPadBottomMargin(0.13);

   TFile fSim("data/attpcsim.root");
   TFile fPRA("data/output_ukf_only.root");
   TFile fIDL("data/output_ukf_truthpra.root", "READ");
   bool haveIdeal = !fIDL.IsZombie() && fIDL.Get("cbmsim") != nullptr;

   auto *tSim = (TTree *)fSim.Get("cbmsim");
   auto *tPRA = (TTree *)fPRA.Get("cbmsim");
   auto *tIDL = haveIdeal ? (TTree *)fIDL.Get("cbmsim") : nullptr;

   TClonesArray *trks = new TClonesArray("AtMCTrack");
   TClonesArray *teP = new TClonesArray("AtTrackingEvent");
   TClonesArray *teI = new TClonesArray("AtTrackingEvent");
   tSim->SetBranchAddress("MCTrack", &trks);
   tPRA->SetBranchAddress("AtTrackingEvent", &teP);
   if (tIDL) tIDL->SetBranchAddress("AtTrackingEvent", &teI);

   TH2F *hTruth = new TH2F("hTruth", ";#theta_{lab} (deg);KE (MeV)", 90, 0., 180., 60, 0., 60.);
   TH2F *hPRA = new TH2F("hPRA", ";#theta_{lab} (deg);KE (MeV)", 90, 0., 180., 60, 0., 60.);
   TH2F *hIdeal = new TH2F("hIdeal", ";#theta_{lab} (deg);KE (MeV)", 90, 0., 180., 60, 0., 60.);

   auto fillFromUkf = [&](TTree *t, TClonesArray *te, TH2F *h) -> int {
      Long64_t n = std::min(tSim->GetEntries(), t->GetEntries());
      int nFit = 0;
      for (Long64_t i = 0; i < n; ++i) {
         tSim->GetEntry(i);
         t->GetEntry(i);
         if (trks->GetEntries() == 0) continue;
         auto *mc = (AtMCTrack *)trks->At(0);
         if (std::abs(mc->GetPdgCode()) != 211) continue;
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
         double thFit = best->GetKinematics().theta * 180. / M_PI;
         if (KEfit < 0 || KEfit > 80) continue;
         h->Fill(thFit, KEfit);
         ++nFit;
      }
      return nFit;
   };

   int nTruth = 0;
   Long64_t n = tSim->GetEntries();
   for (Long64_t i = 0; i < n; ++i) {
      tSim->GetEntry(i);
      if (trks->GetEntries() == 0) continue;
      auto *mc = (AtMCTrack *)trks->At(0);
      if (std::abs(mc->GetPdgCode()) != 211) continue;
      double KEmc = (mc->GetEnergy() - mc->GetMass()) * 1000.;
      double pmc = std::sqrt(mc->GetPx() * mc->GetPx() + mc->GetPy() * mc->GetPy() + mc->GetPz() * mc->GetPz());
      double thMC = std::acos(mc->GetPz() / pmc) * 180. / M_PI;
      hTruth->Fill(thMC, KEmc);
      ++nTruth;
   }
   int nPRA = fillFromUkf(tPRA, teP, hPRA);
   int nIdeal = haveIdeal ? fillFromUkf(tIDL, teI, hIdeal) : 0;

   // Side-by-side: truth | PRA reco | (truth-PRA reco)
   int nPanels = haveIdeal ? 3 : 2;
   auto *c = new TCanvas("c", "pi-TPC kinematics", 600 * nPanels, 600);
   c->Divide(nPanels, 1, 0.005, 0.02);

   c->cd(1);
   gPad->SetGrid();
   hTruth->SetTitle(Form("MC truth (%d);#theta_{lab} (deg);KE (MeV)", nTruth));
   hTruth->SetMarkerStyle(20);
   hTruth->SetMarkerSize(0.4);
   hTruth->SetMarkerColor(kBlack);
   hTruth->Draw();

   c->cd(2);
   gPad->SetGrid();
   hPRA->SetTitle(Form("UKF (PRA seed) (%d);#theta_{lab} (deg);KE (MeV)", nPRA));
   hPRA->SetMarkerStyle(20);
   hPRA->SetMarkerSize(0.4);
   hPRA->SetMarkerColor(kAzure + 2);
   hPRA->Draw();

   if (haveIdeal) {
      c->cd(3);
      gPad->SetGrid();
      hIdeal->SetTitle(Form("UKF (truth seed = intrinsic ceiling) (%d);#theta_{lab} (deg);KE (MeV)", nIdeal));
      hIdeal->SetMarkerStyle(20);
      hIdeal->SetMarkerSize(0.4);
      hIdeal->SetMarkerColor(kRed + 1);
      hIdeal->Draw();
   }

   c->SaveAs("data/kinematics.png");
   c->SaveAs("data/kinematics.pdf");

   // Overlay: truth (gray fill) + PRA (blue) + ideal (red), single panel
   auto *cOv = new TCanvas("cOv", "pi-TPC kinematics overlay", 1000, 700);
   gPad->SetGrid();
   hTruth->SetTitle("Truth (black) vs UKF [PRA seed] (blue) vs UKF [truth seed] (red);#theta_{lab} (deg);KE (MeV)");
   hTruth->SetMarkerSize(0.5);
   hTruth->SetMarkerColor(kGray + 2);
   hTruth->Draw();
   hPRA->SetMarkerSize(0.4);
   hPRA->SetMarkerColor(kAzure + 2);
   hPRA->Draw("same");
   if (haveIdeal) {
      hIdeal->SetMarkerSize(0.4);
      hIdeal->SetMarkerColor(kRed + 1);
      hIdeal->Draw("same");
   }
   auto *leg = new TLegend(0.55, 0.74, 0.88, 0.88);
   leg->SetBorderSize(0);
   leg->SetFillStyle(0);
   leg->AddEntry(hTruth, Form("MC truth (%d)", nTruth), "p");
   leg->AddEntry(hPRA, Form("UKF PRA seed (%d)", nPRA), "p");
   if (haveIdeal) leg->AddEntry(hIdeal, Form("UKF truth seed (%d)", nIdeal), "p");
   leg->Draw();
   cOv->SaveAs("data/kinematics_overlay.png");
   cOv->SaveAs("data/kinematics_overlay.pdf");

   std::cout << "Wrote data/kinematics.{png,pdf} and kinematics_overlay.{png,pdf}\n";
   std::cout << "MC truth pions: " << nTruth << "\n";
   std::cout << "Reco PRA: " << nPRA << " (" << 100. * nPRA / std::max(1, nTruth) << "%)\n";
   if (haveIdeal) std::cout << "Reco truth seed: " << nIdeal << " (" << 100. * nIdeal / std::max(1, nTruth) << "%)\n";
}
