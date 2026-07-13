/// @file make_3d_display.C
/// @brief 3D event display of reconstructed PUMA tracks: hits in (x, y, z_drift),
///        coloured per pattern-recognition track, inside a wireframe of the annular
///        TPC volume (inner R=62.9, outer R=121.1 mm; z = drift depth). Shows how
///        the pad (x,y) + drift-time (z) give the full 3D track.
/// Run: root -b -q 'make_3d_display.C("data/reco_pid_base.root",8)'
static void ring(double R, double z, Color_t col, int style = 1)
{
   auto *l = new TPolyLine3D(65);
   for (int i = 0; i <= 64; ++i) { double a = 2 * M_PI * i / 64; l->SetPoint(i, R * cos(a), R * sin(a), z); }
   l->SetLineColor(col); l->SetLineStyle(style); l->Draw();
}

void make_3d_display(TString file = "data/reco_pid_base.root", Long64_t iev = 8,
                     TString out = "/Users/quantumlab/fair_install/puma_slides/figs/track_3d.png")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62, "xyz"); gStyle->SetTitleFont(62, "xyz");

   TFile f(file); auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pat = nullptr; t->SetBranchAddress("AtPatternEvent", &pat);
   t->GetEntry(iev);
   auto *pe = (AtPatternEvent *)(pat ? pat->At(0) : nullptr);
   if (!pe) { printf("no pattern event\n"); return; }
   auto &tracks = pe->GetTrackCand();

   // z range of the hits
   double zmin = 1e9, zmax = -1e9;
   for (auto &tr : tracks) for (auto &h : tr.GetHitArray()) { double z = h->GetPosition().Z(); zmin = std::min(zmin, z); zmax = std::max(zmax, z); }
   if (zmin > zmax) { zmin = 0; zmax = 150; }
   double zlo = std::max(0.0, zmin - 15), zhi = zmax + 15;

   auto *c = new TCanvas("c3d", "", 1000, 900);
   auto *fr = new TH3F("fr", "PUMA 3D event display;x [mm];y [mm];z_{drift} [mm]",
                       1, -130, 130, 1, -130, 130, 1, zlo, zhi);
   fr->Draw();

   // TPC volume wireframe: outer + inner annulus at the two z ends + vertical struts
   for (double z : {zlo, zhi}) { ring(121.1, z, kGray + 1); ring(62.9, z, kGray + 1); }
   for (int i = 0; i < 8; ++i) { double a = 2 * M_PI * i / 8;
      for (double R : {62.9, 121.1}) { auto *v = new TPolyLine3D(2);
         v->SetPoint(0, R * cos(a), R * sin(a), zlo); v->SetPoint(1, R * cos(a), R * sin(a), zhi);
         v->SetLineColor(kGray); v->SetLineStyle(3); v->Draw(); } }
   // a few pad-plane rings at the pad plane (z=zlo) for scale
   double r2s = (121.1 * 121.1 - 62.9 * 62.9) / 16, rr = 62.9;
   for (int i = 0; i < 16; ++i) { rr = std::sqrt(rr * rr + r2s); ring(rr, zlo, kGray, 1); }

   // tracks: one colour per PRA track candidate
   Color_t cols[] = {kRed + 1, kAzure + 2, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kViolet};
   int it = 0;
   for (auto &tr : tracks) {
      auto *pm = new TPolyMarker3D(tr.GetHitArray().size());
      int n = 0; for (auto &h : tr.GetHitArray()) { const auto &p = h->GetPosition(); pm->SetPoint(n++, p.X(), p.Y(), p.Z()); }
      pm->SetMarkerStyle(20); pm->SetMarkerSize(1.0); pm->SetMarkerColor(cols[it % 7]); pm->Draw();
      it++;
   }
   gPad->SetTheta(28); gPad->SetPhi(35);
   auto *tx = new TLatex(); tx->SetNDC(); tx->SetTextFont(62); tx->SetTextSize(0.03);
   tx->DrawLatex(0.13, 0.90, Form("event %lld  -  %d reconstructed tracks", iev, (int)tracks.size()));
   tx->DrawLatex(0.13, 0.86, "colour = track;  z = drift depth (from pulse time)");
   c->SaveAs(out);
   printf("TRACK3D_DONE event %lld, %d tracks, z[%.0f,%.0f]\n", iev, (int)tracks.size(), zlo, zhi);
}
