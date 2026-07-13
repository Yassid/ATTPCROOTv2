/// @file make_gallery.C
/// @brief Gallery of cleanly reconstructed PUMA events: a grid of pad-plane views,
///        each fired pad drawn as its real annular sector coloured by the PRA track
///        it belongs to. Selects "good" events (2-5 tracks, each with enough hits).
/// Run: root -b -q 'make_gallery.C("data/reco_pid_base.root",9)'
void make_gallery(TString file = "data/reco_pid_base.root", int nShow = 9, int minHitsTrk = 12,
                  TString out = "/Users/quantumlab/fair_install/puma_slides/figs/gallery.png")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62, "xyz"); gStyle->SetTitleFont(62, "xyz");
   gStyle->SetTitleSize(0.07, "t"); gStyle->SetTitleFontSize(0.07);

   // pad geometry (AtTpcPUMAMap)
   const double rIn = 62.9, rOut = 121.1; const int nR = 16, nP = 256;
   std::vector<double> rr(nR + 1); rr[0] = rIn; double r2s = (rOut * rOut - rIn * rIn) / nR;
   for (int i = 1; i <= nR; ++i) rr[i] = std::sqrt(rr[i - 1] * rr[i - 1] + r2s);
   const double dPhi = 2 * M_PI / nP;
   auto padPoly = [&](int pad, Color_t col) {
      int ring = pad / nP, sec = pad % nP; if (ring < 0 || ring >= nR) return;
      double ri = rr[ring], ro = rr[ring + 1], p0 = (sec - 0.5) * dPhi, p1 = (sec + 0.5) * dPhi;
      auto *g = new TGraph(); int k = 0;
      for (int i = 0; i <= 3; ++i) { double p = p0 + (p1 - p0) * i / 3; g->SetPoint(k++, ri * cos(p), ri * sin(p)); }
      for (int i = 0; i <= 3; ++i) { double p = p1 + (p0 - p1) * i / 3; g->SetPoint(k++, ro * cos(p), ro * sin(p)); }
      g->SetPoint(k, ri * cos(p0), ri * sin(p0));
      g->SetFillColor(col); g->SetLineColor(col); g->Draw("f same");
   };

   TFile f(file); auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pat = nullptr; t->SetBranchAddress("AtPatternEvent", &pat);

   // select good events
   std::vector<Long64_t> good;
   for (Long64_t i = 0; i < t->GetEntries() && (int)good.size() < nShow; ++i) {
      t->GetEntry(i);
      auto *pe = (AtPatternEvent *)(pat ? pat->At(0) : nullptr); if (!pe) continue;
      auto &tc = pe->GetTrackCand(); int nt = tc.size();
      if (nt < 2 || nt > 5) continue;
      bool ok = true; for (auto &tr : tc) if ((int)tr.GetHitArray().size() < minHitsTrk) { ok = false; break; }
      if (ok) good.push_back(i);
   }
   int ng = good.size(); int nc = (ng <= 4) ? 2 : 3; int nrpan = (ng + nc - 1) / nc;
   auto *c = new TCanvas("cg", "", 360 * nc, 360 * nrpan);
   c->Divide(nc, nrpan, 0.002, 0.012);
   Color_t cols[] = {kRed + 1, kAzure + 2, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet};
   for (int k = 0; k < ng; ++k) {
      c->cd(k + 1); gPad->SetMargin(0.03, 0.03, 0.03, 0.03); // equal margins -> square frame -> round rings
      t->GetEntry(good[k]);
      auto *pe = (AtPatternEvent *)pat->At(0); auto &tc = pe->GetTrackCand();
      auto *fr = new TH2F(Form("fr%d", k), "", 1, -136, 136, 1, -136, 136); // equal ranges -> round rings
      fr->GetXaxis()->SetLabelSize(0); fr->GetYaxis()->SetLabelSize(0);
      fr->GetXaxis()->SetTickLength(0); fr->GetYaxis()->SetTickLength(0); fr->Draw();
      auto *o = new TEllipse(0, 0, rOut, rOut); o->SetFillStyle(0); o->SetLineColor(16); o->Draw();
      auto *in = new TEllipse(0, 0, rIn, rIn); in->SetFillStyle(0); in->SetLineColor(16); in->Draw();
      int it = 0; for (auto &tr : tc) { Color_t col = cols[it % 7];
         for (auto &h : tr.GetHitArray()) padPoly(h->GetPadNum(), col); it++; }
      auto *lab = new TLatex(0, 127, Form("event %lld  #upoint  %d tracks", good[k], (int)tc.size()));
      lab->SetTextFont(62); lab->SetTextSize(0.052); lab->SetTextAlign(22); lab->Draw();
   }
   c->SaveAs(out);
   printf("GALLERY_DONE %d events\n", ng);
}
