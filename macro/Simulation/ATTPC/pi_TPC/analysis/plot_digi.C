/// @file plot_digi.C
/// @brief Visualise the digi output: pad-plane occupancy + per-event displays.
///
/// Reads ./data/output_digi.root (AtEventH = PSA hits, AtPatternEvent = PRA tracks)
/// and ./data/attpcsim.root (MC truth points). Produces:
///   data/digi_overview.png  -- pad-plane hit map (all events, log scale) + Z range
///   data/digi_events.png    -- 4-event x,y vs z displays with MC truth overlay
///
/// Run: root -b -q plot_digi.C

void plot_digi(int eventDisplayCount = 4)
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kViridis);

   TFile fDigi("data/output_digi.root");
   TFile fSim("data/attpcsim.root");
   auto *tDigi = (TTree *)fDigi.Get("cbmsim");
   auto *tSim = (TTree *)fSim.Get("cbmsim");
   if (!tDigi || !tSim) {
      std::cout << "missing trees\n";
      return;
   }

   TClonesArray *evt = nullptr;
   TClonesArray *peArr = nullptr;
   TClonesArray *mcPts = new TClonesArray("AtMCPoint");
   tDigi->SetBranchAddress("AtEventH", &evt);
   if (tDigi->GetBranch("AtPatternEvent")) tDigi->SetBranchAddress("AtPatternEvent", &peArr);
   tSim->SetBranchAddress("AtTpcPoint", &mcPts);

   // ----- Aggregate occupancy --------------------------------------------
   TH2F *hOcc = new TH2F("hOcc", "Pad-plane hit density (all events);x (mm);y (mm)", 110, -110., 110., 110, -110., 110.);
   TH2F *hZxy = new TH2F("hZxy", "Mean z (mm) per pad;x (mm);y (mm);<z> (mm)", 110, -110., 110., 110, -110., 110.);
   TH2F *hCnt = new TH2F("hCnt", "", 110, -110., 110., 110, -110., 110.);
   TH1F *hZ = new TH1F("hZ", "Hit z distribution;z (mm);counts", 100, 0., 1000.);
   TH1F *hQ = new TH1F("hQ", "Hit charge;charge (a.u.);counts", 100, 0., 5000.);

   Long64_t n = std::min(tDigi->GetEntries(), tSim->GetEntries());
   for (Long64_t i = 0; i < n; ++i) {
      tDigi->GetEntry(i);
      auto *e = (AtEvent *)evt->At(0);
      if (!e) continue;
      for (const auto &h : e->GetHits()) {
         const auto &p = h->GetPosition();
         hOcc->Fill(p.X(), p.Y());
         hZxy->Fill(p.X(), p.Y(), p.Z());
         hCnt->Fill(p.X(), p.Y());
         hZ->Fill(p.Z());
         hQ->Fill(h->GetCharge());
      }
   }
   hZxy->Divide(hCnt);

   auto *cAgg = new TCanvas("cAgg", "Digi overview", 1500, 1000);
   cAgg->Divide(2, 2);
   cAgg->cd(1);
   gPad->SetLogz();
   gPad->SetGrid();
   hOcc->Draw("colz");
   cAgg->cd(2);
   gPad->SetGrid();
   hZxy->Draw("colz");
   cAgg->cd(3);
   gPad->SetGrid();
   hZ->SetFillColorAlpha(kAzure - 4, 0.6);
   hZ->SetLineColor(kAzure + 2);
   hZ->Draw("hist");
   cAgg->cd(4);
   gPad->SetGrid();
   gPad->SetLogy();
   hQ->SetFillColorAlpha(kSpring - 7, 0.6);
   hQ->SetLineColor(kGreen + 2);
   hQ->Draw("hist");
   cAgg->SaveAs("data/digi_overview.png");
   cAgg->SaveAs("data/digi_overview.pdf");

   // ----- Per-event displays --------------------------------------------
   auto *cEv = new TCanvas("cEv", "Event displays", 1600, 1200);
   cEv->Divide(2, eventDisplayCount / 2 + (eventDisplayCount % 2));

   int picked = 0;
   for (Long64_t i = 0; i < n && picked < eventDisplayCount; ++i) {
      tDigi->GetEntry(i);
      tSim->GetEntry(i);
      auto *e = (AtEvent *)evt->At(0);
      if (!e || e->GetHits().empty()) continue;
      picked++;

      // Hits as graph in (x, y) and (z, x).
      TGraph *gHitsXY = new TGraph();
      TGraph *gHitsZX = new TGraph();
      int idx = 0;
      for (const auto &h : e->GetHits()) {
         const auto &p = h->GetPosition();
         gHitsXY->SetPoint(idx, p.X(), p.Y());
         gHitsZX->SetPoint(idx, p.Z(), std::sqrt(p.X() * p.X() + p.Y() * p.Y()) * (p.X() < 0 ? -1 : 1));
         idx++;
      }

      TGraph *gMCXY = new TGraph();
      TGraph *gMCZX = new TGraph();
      int mIdx = 0;
      for (int k = 0; k < mcPts->GetEntries(); ++k) {
         auto *mp = (AtMCPoint *)mcPts->At(k);
         double mx = mp->GetX() * 10., my = mp->GetY() * 10., mz = mp->GetZ() * 10.;
         gMCXY->SetPoint(mIdx, mx, my);
         gMCZX->SetPoint(mIdx, mz, std::sqrt(mx * mx + my * my) * (mx < 0 ? -1 : 1));
         mIdx++;
      }

      cEv->cd(picked);
      gPad->SetGrid();
      gHitsXY->SetMarkerStyle(20);
      gHitsXY->SetMarkerSize(0.5);
      gHitsXY->SetMarkerColor(kAzure + 2);
      gHitsXY->SetTitle(Form("Event %lld;x (mm);y (mm)", i));
      gHitsXY->Draw("AP");
      // Set axis range
      gHitsXY->GetXaxis()->SetLimits(-110., 110.);
      gHitsXY->SetMinimum(-110.);
      gHitsXY->SetMaximum(110.);
      gMCXY->SetMarkerStyle(2);
      gMCXY->SetMarkerSize(0.4);
      gMCXY->SetMarkerColor(kRed);
      gMCXY->Draw("P same");
      // Truth vertex marker
      TMarker *vtxMk = new TMarker(0., 0., 29);
      vtxMk->SetMarkerColor(kRed + 1);
      vtxMk->SetMarkerSize(2.0);
      vtxMk->Draw("same");
   }
   cEv->SaveAs("data/digi_events.png");
   cEv->SaveAs("data/digi_events.pdf");

   std::cout << "Wrote data/digi_overview.{png,pdf} and data/digi_events.{png,pdf}\n";
}
