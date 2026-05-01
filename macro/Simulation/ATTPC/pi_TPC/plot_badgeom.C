/// @file plot_badgeom.C
/// @brief Show events that have a PRA TrackCand but no UKF fit — i.e. the
/// 19 "bad-geom" drops where AtFitterUKF rejected on NaN theta or
/// GeoRadius <= 0. Plots MC truth + digi hits + PRA-circle annotation.
/// Run: root -b -q plot_badgeom.C

void plot_badgeom(int nShow = 6)
{
   gStyle->SetOptStat(0);

   TFile fSim("data/attpcsim.root");
   TFile fDigi("data/output_digi.root");
   TFile fUKF("data/output_ukf_only.root");
   auto *tSim = (TTree *)fSim.Get("cbmsim");
   auto *tDigi = (TTree *)fDigi.Get("cbmsim");
   auto *tUKF = (TTree *)fUKF.Get("cbmsim");

   TClonesArray *mcPts = new TClonesArray("AtMCPoint");
   TClonesArray *mcTracks = new TClonesArray("AtMCTrack");
   TClonesArray *evt = nullptr;
   TClonesArray *peArr = nullptr;
   TClonesArray *teArr = nullptr;
   tSim->SetBranchAddress("AtTpcPoint", &mcPts);
   tSim->SetBranchAddress("MCTrack", &mcTracks);
   tDigi->SetBranchAddress("AtEventH", &evt);
   tDigi->SetBranchAddress("AtPatternEvent", &peArr);
   tUKF->SetBranchAddress("AtTrackingEvent", &teArr);

   auto *c = new TCanvas("c", "Bad-geom events", 1500, 1000);
   int rows = (nShow + 2) / 3;
   c->Divide(3, rows);

   Long64_t n = std::min({tSim->GetEntries(), tDigi->GetEntries(), tUKF->GetEntries()});
   int picked = 0;

   for (Long64_t i = 0; i < n && picked < nShow; ++i) {
      tDigi->GetEntry(i);
      tUKF->GetEntry(i);
      if (peArr->GetEntries() == 0 || teArr->GetEntries() == 0) continue;
      auto *pe = (AtPatternEvent *)peArr->At(0);
      auto *te = (AtTrackingEvent *)teArr->At(0);
      if (pe->GetTrackCand().empty() || !te->GetFittedTracks().empty()) continue;
      // Pattern exists but UKF rejected. The remaining drops should all be bad-geom.

      tSim->GetEntry(i);
      auto *e = (AtEvent *)evt->At(0);
      auto *mc = (AtMCTrack *)mcTracks->At(0);
      if (!mc) continue;

      double KEmc = (mc->GetEnergy() - mc->GetMass()) * 1000.;

      // PRA circle / track diagnostics
      auto &cand = pe->GetTrackCand();
      AtTrack &tr = cand.front();
      double praTheta = tr.GetGeoTheta() * 180. / M_PI;
      double praRadius = tr.GetGeoRadius();
      auto praCenter = tr.GetGeoCenter();
      bool nanTheta = std::isnan(praTheta) || std::isnan(tr.GetGeoTheta());
      bool badRadius = std::isnan(praRadius) || praRadius <= 0;

      // MC truth
      TGraph *gMC = new TGraph();
      int mIdx = 0;
      for (int k = 0; k < mcPts->GetEntries(); ++k) {
         auto *mp = (AtMCPoint *)mcPts->At(k);
         gMC->SetPoint(mIdx++, mp->GetX() * 10., mp->GetY() * 10.);
      }
      // Digi hits
      TGraph *gHit = new TGraph();
      int hIdx = 0;
      for (const auto &h : e->GetHits()) {
         const auto &p = h->GetPosition();
         gHit->SetPoint(hIdx++, p.X(), p.Y());
      }

      picked++;
      c->cd(picked);
      gPad->SetGrid();
      gMC->SetMarkerStyle(7);
      gMC->SetMarkerColor(kRed);
      gMC->SetTitle(Form("Event %lld  KE=%.1f MeV  PRA: #theta=%.1f, R=%.1f mm  %s%s;x (mm);y (mm)", i, KEmc,
                         praTheta, praRadius, nanTheta ? "[NaN-theta]" : "", badRadius ? "[bad-R]" : ""));
      gMC->Draw("AP");
      gMC->GetXaxis()->SetLimits(-150., 150.);
      gMC->SetMinimum(-150.);
      gMC->SetMaximum(150.);
      if (hIdx > 0) {
         gHit->SetMarkerStyle(20);
         gHit->SetMarkerSize(0.5);
         gHit->SetMarkerColor(kAzure + 2);
         gHit->Draw("P same");
      }
      // PRA circle (if valid)
      if (!badRadius && std::isfinite(praCenter.first) && std::isfinite(praCenter.second)) {
         TEllipse *circ = new TEllipse(praCenter.first, praCenter.second, praRadius);
         circ->SetFillStyle(0);
         circ->SetLineColor(kMagenta + 2);
         circ->SetLineWidth(2);
         circ->Draw("same");
      }
      // Pad-plane outline
      TBox *padBox = new TBox(-100., -100., 100., 100.);
      padBox->SetFillStyle(0);
      padBox->SetLineColor(kGreen + 2);
      padBox->SetLineWidth(2);
      padBox->Draw("same");
      // Vertex marker
      TMarker *vtxMk = new TMarker(0., 0., 29);
      vtxMk->SetMarkerColor(kRed + 1);
      vtxMk->SetMarkerSize(2.0);
      vtxMk->Draw("same");
   }

   c->SaveAs("data/digi_badgeom.png");
   c->SaveAs("data/digi_badgeom.pdf");
   std::cout << "Wrote data/digi_badgeom.{png,pdf}  (showed " << picked << " bad-geom events)\n";
}
