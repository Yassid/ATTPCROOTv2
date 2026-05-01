/// @file plot_dropped.C
/// @brief Show events where the UKF found NO fitted track. Plots MC truth
/// trajectory in (x,y) and digi hits with the 200x200 mm² pad-plane outline.
/// Run: root -b -q plot_dropped.C

void plot_dropped(int nShow = 6)
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
   TClonesArray *teArr = nullptr;
   tSim->SetBranchAddress("AtTpcPoint", &mcPts);
   tSim->SetBranchAddress("MCTrack", &mcTracks);
   tDigi->SetBranchAddress("AtEventH", &evt);
   tUKF->SetBranchAddress("AtTrackingEvent", &teArr);

   auto *c = new TCanvas("c", "Dropped events", 1500, 1000);
   int rows = (nShow + 2) / 3;
   c->Divide(3, rows);

   Long64_t n = std::min({tSim->GetEntries(), tDigi->GetEntries(), tUKF->GetEntries()});
   int picked = 0;

   for (Long64_t i = 0; i < n && picked < nShow; ++i) {
      tUKF->GetEntry(i);
      if (!teArr || teArr->GetEntries() == 0) continue;
      auto *te = (AtTrackingEvent *)teArr->At(0);
      if (!te) continue;
      // Pick events with NO fitted track
      if (!te->GetFittedTracks().empty()) continue;

      tSim->GetEntry(i);
      tDigi->GetEntry(i);
      auto *e = (AtEvent *)evt->At(0);
      auto *mc = (AtMCTrack *)mcTracks->At(0);
      if (!mc) continue;

      double KEmc = (mc->GetEnergy() - mc->GetMass()) * 1000.;

      // MC truth points
      TGraph *gMC = new TGraph();
      int mIdx = 0;
      for (int k = 0; k < mcPts->GetEntries(); ++k) {
         auto *mp = (AtMCPoint *)mcPts->At(k);
         gMC->SetPoint(mIdx++, mp->GetX() * 10., mp->GetY() * 10.);
      }
      // Digi hits (PSA)
      TGraph *gHit = new TGraph();
      int hIdx = 0;
      if (e) {
         for (const auto &h : e->GetHits()) {
            const auto &p = h->GetPosition();
            gHit->SetPoint(hIdx++, p.X(), p.Y());
         }
      }

      picked++;
      c->cd(picked);
      gPad->SetGrid();

      // Draw MC truth (red line) first, then hits (blue) on top
      gMC->SetMarkerStyle(7);
      gMC->SetMarkerColor(kRed);
      gMC->SetTitle(
         Form("Event %lld  KE_{MC}=%.1f MeV  N_{hit}=%d  N_{MC}=%d;x (mm);y (mm)", i, KEmc, hIdx, mIdx));
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

      // Pad plane outline (200x200 mm² square)
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

   if (picked == 0) {
      std::cout << "No dropped events found in first " << n << " entries.\n";
   }
   c->SaveAs("data/digi_dropped.png");
   c->SaveAs("data/digi_dropped.pdf");
   std::cout << "Wrote data/digi_dropped.{png,pdf}  (showed " << picked << " dropped events)\n";
}
