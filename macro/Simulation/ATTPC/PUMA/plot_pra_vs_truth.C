/// @file plot_pra_vs_truth.C
/// @brief Per-event display: MC truth hits coloured by primary particle vs
/// PRA track candidates coloured by track-id. Picks events with >=2 PRA
/// tracks so we can see how TripletClust separates (or merges/splits) the
/// PUMA multi-particle final state.
///
/// Run: root -b -q plot_pra_vs_truth.C

void plot_pra_vs_truth(int nShow = 6)
{
   gStyle->SetOptStat(0);

   TFile fSim("data/attpcsim.root");
   TFile fDigi("data/output_digi.root");
   auto *tSim = (TTree *)fSim.Get("cbmsim");
   auto *tDigi = (TTree *)fDigi.Get("cbmsim");
   if (!tSim || !tDigi) {
      std::cout << "missing trees\n";
      return;
   }

   TClonesArray *mcPts = new TClonesArray("AtMCPoint");
   TClonesArray *peArr = nullptr;
   tSim->SetBranchAddress("AtTpcPoint", &mcPts);
   tDigi->SetBranchAddress("AtPatternEvent", &peArr);

   const Color_t mcCol[] = {kRed, kBlue, kGreen + 2, kMagenta, kOrange + 7, kCyan + 2, kViolet + 1};
   const Color_t praCol[] = {kBlack, kRed + 1, kBlue + 1, kGreen + 3, kMagenta + 2, kOrange + 8, kCyan + 3};

   auto *cXY = new TCanvas("cXY", "PRA vs MC truth (x,y)", 1500, 1000);
   cXY->Divide(3, 2);

   Long64_t n = std::min(tSim->GetEntries(), tDigi->GetEntries());
   int picked = 0;
   for (Long64_t i = 0; i < n && picked < nShow; ++i) {
      tDigi->GetEntry(i);
      if (peArr->GetEntries() == 0) continue;
      auto *pe = (AtPatternEvent *)peArr->At(0);
      auto &cand = pe->GetTrackCand();
      if (cand.size() < 2) continue;
      tSim->GetEntry(i);

      picked++;
      cXY->cd(picked);
      gPad->SetGrid();

      // Frame
      auto *hFrame = new TH2F(Form("hF%d", picked),
                              Form("Event %lld  N_{PRA}=%d  N_{MC tracks}=%d;x (mm);y (mm)", i, (int)cand.size(),
                                   mcPts->GetEntries()),
                              10, -150., 150., 10, -150., 150.);
      hFrame->Draw();

      // Pad-plane outline
      TBox *padBox = new TBox(-100., -100., 100., 100.);
      padBox->SetFillStyle(0);
      padBox->SetLineColor(kGray + 1);
      padBox->SetLineWidth(1);
      padBox->Draw("same");

      // MC truth points coloured by trackID
      std::map<int, TGraph *> mcByTrk;
      for (int k = 0; k < mcPts->GetEntries(); ++k) {
         auto *mp = (AtMCPoint *)mcPts->At(k);
         int tid = mp->GetTrackID();
         if (mcByTrk.count(tid) == 0) {
            mcByTrk[tid] = new TGraph();
            mcByTrk[tid]->SetMarkerStyle(7);
            mcByTrk[tid]->SetMarkerColor(mcCol[mcByTrk.size() % 7]);
         }
         auto *g = mcByTrk[tid];
         g->SetPoint(g->GetN(), mp->GetX() * 10., mp->GetY() * 10.);
      }
      for (auto &kv : mcByTrk) kv.second->Draw("P same");

      // PRA candidate hits coloured by track index
      int praIdx = 0;
      for (auto &tr : cand) {
         auto &hits = tr.GetHitArray();
         if (hits.empty()) continue;
         auto *g = new TGraph();
         g->SetMarkerStyle(24); // open circle
         g->SetMarkerSize(0.5);
         g->SetMarkerColor(praCol[praIdx % 7]);
         for (auto &hp : hits) {
            const auto &p = hp->GetPosition();
            g->SetPoint(g->GetN(), p.X(), p.Y());
         }
         g->Draw("P same");
         praIdx++;
      }

      // Vertex marker
      TMarker *vtxMk = new TMarker(0., 0., 29);
      vtxMk->SetMarkerColor(kRed + 2);
      vtxMk->SetMarkerSize(1.8);
      vtxMk->Draw("same");
   }

   cXY->SaveAs("data/pra_vs_truth.png");
   cXY->SaveAs("data/pra_vs_truth.pdf");
   std::cout << "Wrote data/pra_vs_truth.{png,pdf}  (" << picked << " events shown)\n";
}
