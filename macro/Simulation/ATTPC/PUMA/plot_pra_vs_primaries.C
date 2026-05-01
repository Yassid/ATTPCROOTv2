/// @file plot_pra_vs_primaries.C
/// @brief Per-event display showing only the PRIMARY K/pi tracks that reach
/// the gas drift volume (>=30 MC points), vs the PRA candidate tracks.
/// This isolates the actual PRA efficiency on reconstructible tracks,
/// without the visual clutter of secondaries / decay products / delta-rays.
///
/// Run: root -b -q plot_pra_vs_primaries.C

void plot_pra_vs_primaries(int nShow = 6, int minHitsPrim = 30)
{
   gStyle->SetOptStat(0);

   TFile fSim("data/attpcsim.root");
   TFile fDigi("data/output_digi.root");
   auto *tSim = (TTree *)fSim.Get("cbmsim");
   auto *tDigi = (TTree *)fDigi.Get("cbmsim");

   TClonesArray *mcPts = new TClonesArray("AtMCPoint");
   TClonesArray *mcTracks = new TClonesArray("AtMCTrack");
   TClonesArray *peArr = nullptr;
   tSim->SetBranchAddress("AtTpcPoint", &mcPts);
   tSim->SetBranchAddress("MCTrack", &mcTracks);
   tDigi->SetBranchAddress("AtPatternEvent", &peArr);

   const Color_t mcCol[] = {kRed, kBlue, kGreen + 2};
   const Color_t praCol[] = {kBlack, kRed + 1, kBlue + 1, kGreen + 3, kMagenta + 2, kOrange + 8};

   auto *cXY = new TCanvas("cXY", "Primaries vs PRA (x,y)", 1500, 1000);
   cXY->Divide(3, 2);

   Long64_t n = std::min(tSim->GetEntries(), tDigi->GetEntries());
   int picked = 0;
   for (Long64_t i = 0; i < n && picked < nShow; ++i) {
      tSim->GetEntry(i);

      // Build per-trackID MC point lists (gas only)
      std::map<int, std::vector<std::pair<double, double>>> primTracks;
      std::map<int, int> hitsByTrack;
      for (int k = 0; k < mcPts->GetEntries(); ++k) {
         auto *mp = (AtMCPoint *)mcPts->At(k);
         hitsByTrack[mp->GetTrackID()]++;
      }

      // Filter to primary K/pi with enough hits
      std::vector<int> keep;
      std::vector<std::string> names;
      for (int k = 0; k < mcTracks->GetEntries(); ++k) {
         auto *tr = (AtMCTrack *)mcTracks->At(k);
         if (tr->GetMotherId() != -1) continue;
         int pdg = tr->GetPdgCode();
         int absPdg = std::abs(pdg);
         if (absPdg != 211 && absPdg != 321) continue;
         if (hitsByTrack[k] < minHitsPrim) continue;
         keep.push_back(k);
         std::string nm = (absPdg == 321) ? "K" : "pi";
         nm += (pdg > 0) ? "+" : "-";
         names.push_back(nm);
      }
      if (keep.empty()) continue;

      tDigi->GetEntry(i);
      // Now collect points for the kept tracks
      for (int kid : keep) primTracks[kid] = {};
      for (int k = 0; k < mcPts->GetEntries(); ++k) {
         auto *mp = (AtMCPoint *)mcPts->At(k);
         int t = mp->GetTrackID();
         if (primTracks.count(t))
            primTracks[t].push_back({mp->GetX() * 10., mp->GetY() * 10.});
      }

      picked++;
      cXY->cd(picked);
      gPad->SetGrid();

      auto *hFrame = new TH2F(Form("hF%d", picked),
                              Form("Event %lld  N_{prim in gas}=%d  N_{PRA}=%d;x (mm);y (mm)", i, (int)keep.size(),
                                   peArr->GetEntries() ? (int)((AtPatternEvent *)peArr->At(0))->GetTrackCand().size()
                                                       : 0),
                              10, -150., 150., 10, -150., 150.);
      hFrame->Draw();

      // Annular drift outline
      TEllipse *rOut = new TEllipse(0, 0, 125.1);
      rOut->SetFillStyle(0);
      rOut->SetLineColor(kGray + 2);
      rOut->Draw("same");
      TEllipse *rIn = new TEllipse(0, 0, 58.5);
      rIn->SetFillStyle(0);
      rIn->SetLineColor(kGray + 2);
      rIn->Draw("same");

      // Pad-plane outline
      TBox *padBox = new TBox(-100., -100., 100., 100.);
      padBox->SetFillStyle(0);
      padBox->SetLineColor(kGray + 1);
      padBox->SetLineStyle(2);
      padBox->Draw("same");

      // Primary tracks (colored by particle, big markers, with label)
      double labelY = 130.;
      int idx = 0;
      for (int kid : keep) {
         auto *g = new TGraph();
         for (auto &p : primTracks[kid]) g->SetPoint(g->GetN(), p.first, p.second);
         g->SetMarkerStyle(20);
         g->SetMarkerSize(0.4);
         g->SetMarkerColor(mcCol[idx % 3]);
         g->Draw("P same");
         TLatex *lab = new TLatex(-145, labelY, Form("%s (%d hits)", names[idx].c_str(), hitsByTrack[kid]));
         lab->SetTextColor(mcCol[idx % 3]);
         lab->SetTextSize(0.035);
         lab->Draw();
         labelY -= 14;
         idx++;
      }

      // PRA candidate tracks (open circles, distinct from primaries)
      if (peArr->GetEntries() > 0) {
         auto *pe = (AtPatternEvent *)peArr->At(0);
         int praIdx = 0;
         for (auto &tr : pe->GetTrackCand()) {
            auto &hits = tr.GetHitArray();
            if (hits.empty()) continue;
            auto *g = new TGraph();
            g->SetMarkerStyle(24);
            g->SetMarkerSize(0.7);
            g->SetMarkerColor(praCol[praIdx % 6]);
            for (auto &hp : hits) {
               const auto &p = hp->GetPosition();
               g->SetPoint(g->GetN(), p.X(), p.Y());
            }
            g->Draw("P same");
            praIdx++;
         }
      }

      // Vertex marker
      TMarker *vtxMk = new TMarker(0., 0., 29);
      vtxMk->SetMarkerColor(kRed + 2);
      vtxMk->SetMarkerSize(1.8);
      vtxMk->Draw("same");
   }

   cXY->SaveAs("data/pra_vs_primaries.png");
   cXY->SaveAs("data/pra_vs_primaries.pdf");
   std::cout << "Wrote data/pra_vs_primaries.{png,pdf}  (" << picked << " events shown)\n";
}
