/// @file make_2d_tracks.C
/// @brief 2D pad-plane view of one event: each fired pad drawn as its real annular
///        sector, coloured by the pattern-recognition track it belongs to. Pairs
///        with make_3d_display.C (same event) to show the 2D<->3D correspondence.
/// Run: root -b -q 'make_2d_tracks.C("data/reco_pid_base.root",8)'
void make_2d_tracks(TString file = "data/reco_pid_base.root", Long64_t iev = 8,
                    TString out = "/Users/quantumlab/fair_install/puma_slides/figs/track_2d.png")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62, "xyz"); gStyle->SetTitleFont(62, "xyz");
   gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

   // pad geometry (mirrors AtTpcPUMAMap): 16 equal-area rings x 256 sectors
   const double rIn = 62.9, rOut = 121.1; const int nR = 16, nP = 256;
   std::vector<double> rr(nR + 1); rr[0] = rIn; double r2s = (rOut * rOut - rIn * rIn) / nR;
   for (int i = 1; i <= nR; ++i) rr[i] = std::sqrt(rr[i - 1] * rr[i - 1] + r2s);
   const double dPhi = 2 * M_PI / nP;
   auto padPoly = [&](int pad, Color_t col) {
      int ring = pad / nP, sec = pad % nP; double ri = rr[ring], ro = rr[ring + 1];
      double p0 = (sec - 0.5) * dPhi, p1 = (sec + 0.5) * dPhi;
      auto *g = new TGraph(); int NA = 3, k = 0;
      for (int i = 0; i <= NA; ++i) { double p = p0 + (p1 - p0) * i / NA; g->SetPoint(k++, ri * cos(p), ri * sin(p)); }
      for (int i = 0; i <= NA; ++i) { double p = p1 + (p0 - p1) * i / NA; g->SetPoint(k++, ro * cos(p), ro * sin(p)); }
      g->SetPoint(k, ri * cos(p0), ri * sin(p0));
      g->SetFillColor(col); g->SetLineColor(col); g->SetLineWidth(1); g->Draw("f same");
   };

   TFile f(file); auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pat = nullptr; t->SetBranchAddress("AtPatternEvent", &pat);
   t->GetEntry(iev);
   auto *pe = (AtPatternEvent *)(pat ? pat->At(0) : nullptr);
   if (!pe) { printf("no pattern event\n"); return; }
   auto &tracks = pe->GetTrackCand();

   auto *c = new TCanvas("c2d", "", 900, 850); c->SetLeftMargin(0.13); c->SetBottomMargin(0.12);
   auto *fr = new TH2F("fr", Form("PUMA pad plane - event %lld  (%d tracks);x [mm];y [mm]", iev, (int)tracks.size()),
                       10, -128, 128, 10, -128, 128);
   fr->Draw();
   // faint ring outlines for context
   for (double r : rr) { auto *e = new TEllipse(0, 0, r, r); e->SetFillStyle(0); e->SetLineColor(17); e->Draw(); }

   Color_t cols[] = {kRed + 1, kAzure + 2, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet};
   int it = 0;
   for (auto &tr : tracks) { Color_t col = cols[it % 7];
      for (auto &h : tr.GetHitArray()) { int pad = h->GetPadNum(); if (pad >= 0 && pad < nR * nP) padPoly(pad, col); }
      it++;
   }
   gPad->RedrawAxis();
   c->SaveAs(out);
   printf("TRACK2D_DONE event %lld, %d tracks\n", iev, (int)tracks.size());
}
