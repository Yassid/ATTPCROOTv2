/// @file pra_examples.C
/// @brief Gallery of pattern-recognition (PRA) examples: per event, the hits of
///        each track candidate (coloured per track) with the PRA circle fit
///        (GetGeoCenter/GetGeoRadius) overlaid. Shows what pattern recognition
///        produces before the Kalman fit.
/// Run: root -b -q 'pra_examples.C("data/hs_reco375.root",12)'
void pra_examples(TString file = "data/hs_reco375.root", int nShow = 12, int minHitsTrk = 10,
                  TString out = "/Users/quantumlab/fair_install/puma_slides/figs/pra_examples.png")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62, "xyz"); gStyle->SetTitleFont(62, "xyz");
   TFile f(file); auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pat = nullptr; t->SetBranchAddress("AtPatternEvent", &pat);
   std::vector<Long64_t> good;
   for (Long64_t i = 0; i < t->GetEntries() && (int)good.size() < nShow; ++i) { t->GetEntry(i);
      auto *pe = (AtPatternEvent *)(pat ? pat->At(0) : nullptr); if (!pe) continue;
      auto &tc = pe->GetTrackCand(); if (tc.size() < 2 || tc.size() > 4) continue;
      bool ok = true; for (auto &tr : tc) if ((int)tr.GetHitArray().size() < minHitsTrk) { ok = false; break; }
      if (ok) good.push_back(i); }
   int ng = good.size(), nc = 4, nr = (ng + nc - 1) / nc;
   Color_t cols[] = {kRed + 1, kAzure + 2, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet};
   auto *c = new TCanvas("cpra", "", 330 * nc, 330 * nr); c->Divide(nc, nr, 0.001, 0.012);
   for (int k = 0; k < ng; ++k) {
      c->cd(k + 1); gPad->SetMargin(0.03, 0.03, 0.03, 0.03); gPad->SetFixedAspectRatio(kTRUE);
      t->GetEntry(good[k]); auto *pe = (AtPatternEvent *)pat->At(0); auto &tc = pe->GetTrackCand();
      auto *fr = new TH2F(Form("f%d", k), "", 1, -136, 136, 1, -136, 136);
      fr->GetXaxis()->SetLabelSize(0); fr->GetYaxis()->SetLabelSize(0);
      fr->GetXaxis()->SetTickLength(0); fr->GetYaxis()->SetTickLength(0); fr->Draw();
      for (double R : {62.9, 121.1}) { auto *e = new TEllipse(0, 0, R, R); e->SetFillStyle(0); e->SetLineColor(17); e->Draw(); }
      int it = 0; for (auto &tr : tc) { Color_t col = cols[it % 7];
         auto *g = new TGraph(); for (auto &h : tr.GetHitArray()) { const auto &p = h->GetPosition(); g->SetPoint(g->GetN(), p.X(), p.Y()); }
         g->SetMarkerStyle(20); g->SetMarkerSize(0.45); g->SetMarkerColor(col); g->Draw("P same");
         // PRA circle fit
         double R = tr.GetGeoRadius(); auto ctr = tr.GetGeoCenter();
         if (R > 20 && R < 4000) { auto *e = new TEllipse(ctr.first, ctr.second, R, R); e->SetFillStyle(0);
            e->SetLineColor(col); e->SetLineStyle(2); e->SetLineWidth(1); e->Draw(); }
         it++; }
      auto *lab = new TLatex(0, 127, Form("event %lld  #upoint  %d tracks", good[k], (int)tc.size()));
      lab->SetTextFont(62); lab->SetTextSize(0.056); lab->SetTextAlign(22); lab->Draw();
   }
   c->SaveAs(out);
   printf("PRA_EXAMPLES_DONE %d events\n", ng);
}
