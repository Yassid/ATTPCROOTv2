/// @file make_padplane.C
/// @brief Draw the PUMA annular pad plane with the ACTUAL pads (16 equal-area
///        rings x 256 azimuthal sectors = 4096 pads) as a TH2Poly, each fired
///        pad filled by its deposited charge (ADC) on a Z colour scale. This
///        replaces the "staggered points" displays so the segmentation is clear.
///        Geometry mirrors AtTpcPUMAMap exactly (R = 62.9-121.1 mm).
/// Run: root -b -q 'make_padplane.C'      // makes padplane_event + padplane_dlc

// Build a TH2Poly whose bin (pad+1) is the annular-sector pad `pad = ring*256+sec`.
TH2Poly *buildPUMApads(double rIn = 62.9, double rOut = 121.1, int nR = 16, int nP = 256)
{
   auto *hp = new TH2Poly();
   hp->SetName("pumaPads");
   double r2s = (rOut * rOut - rIn * rIn) / nR;
   std::vector<double> rr(nR + 1);
   rr[0] = rIn;
   for (int i = 1; i <= nR; ++i)
      rr[i] = std::sqrt(rr[i - 1] * rr[i - 1] + r2s);
   double dPhi = 2. * M_PI / nP;
   const int NA = 3; // arc subdivisions per pad edge (smooth curved sectors)
   for (int ring = 0; ring < nR; ++ring) {
      for (int sec = 0; sec < nP; ++sec) {
         double ri = rr[ring], ro = rr[ring + 1];
         double p0 = (sec - 0.5) * dPhi, p1 = (sec + 0.5) * dPhi;
         std::vector<double> x, y;
         for (int k = 0; k <= NA; ++k) { double p = p0 + (p1 - p0) * k / NA; x.push_back(ri * cos(p)); y.push_back(ri * sin(p)); }
         for (int k = 0; k <= NA; ++k) { double p = p1 + (p0 - p1) * k / NA; x.push_back(ro * cos(p)); y.push_back(ro * sin(p)); }
         x.push_back(x[0]); y.push_back(y[0]);
         hp->AddBin(x.size(), x.data(), y.data());
      }
   }
   return hp;
}

// Fill the pad plane from one AtEvent (max charge per pad). Returns the (x,y) of
// the single highest-charge pad (for zoom centring).
static void fillEvent(TTree *t, TClonesArray *ev, Long64_t ie, TH2Poly *hp, int &nh, double &cx, double &cy)
{
   hp->ClearBinContents();
   t->GetEntry(ie);
   auto *e = (AtEvent *)(ev ? ev->At(0) : nullptr);
   nh = e ? e->GetNumHits() : 0;
   double qmax = -1;
   cx = 0; cy = 0;
   for (int j = 0; j < nh; ++j) {
      auto &h = e->GetHit(j);
      int pad = h.GetPadNum();
      double q = h.GetCharge();
      if (pad >= 0 && pad < hp->GetNumberOfBins() && q > hp->GetBinContent(pad + 1))
         hp->SetBinContent(pad + 1, q);
      if (q > qmax) { qmax = q; cx = h.GetPosition().X(); cy = h.GetPosition().Y(); }
   }
}

static void styleZ(TH2Poly *hp, const char *title)
{
   hp->SetTitle(Form("%s;x [mm];y [mm]", title));
   hp->GetZaxis()->SetTitle("pad charge [a.u.]");
   hp->GetZaxis()->SetTitleOffset(1.3);
   hp->GetZaxis()->SetLabelSize(0.028);
   hp->SetMinimum(1e-6);
}

// Full plane: only fired pads filled (COLZ, no dense grid) + 16 ring circles for structure.
static void drawFull(TH2Poly *hp, int nh, const std::vector<double> &rr, const char *title)
{
   styleZ(hp, Form("%s  (%d pads)", title, nh));
   hp->GetXaxis()->SetLimits(-128, 128);
   hp->GetYaxis()->SetRangeUser(-128, 128);
   hp->SetLineColor(kGray + 2);
   hp->Draw("COL Z");
   for (double r : rr) { auto *e = new TEllipse(0, 0, r, r); e->SetFillStyle(0); e->SetLineColor(17); e->Draw(); }
}

// Zoom: a +/-win window around (cx,cy) with COLZ L -> individual pad sectors resolved.
static void drawZoom(TH2Poly *hp, int nh, double cx, double cy, double win, const char *title)
{
   styleZ(hp, Form("%s  (%d pads)", title, nh));
   hp->GetXaxis()->SetLimits(cx - win, cx + win);
   hp->GetYaxis()->SetRangeUser(cy - win, cy + win);
   hp->SetLineColor(15);
   hp->SetLineWidth(1);
   hp->Draw("COLZ L");
}

void make_padplane(TString baseFile = "data/reco_pi_base.root", TString dlcFile = "data/reco_pi_dlc.root",
                   TString outdir = "/Users/quantumlab/fair_install/puma_slides/figs")
{
   gStyle->SetOptStat(0);
   gStyle->SetTextFont(62);
   gStyle->SetLabelFont(62, "xyz");
   gStyle->SetTitleFont(62, "xyz");
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(99);

   // ring radii for the structure circles
   double rIn = 62.9, rOut = 121.1; int nR = 16;
   std::vector<double> rr(nR + 1);
   rr[0] = rIn; double r2s = (rOut * rOut - rIn * rIn) / nR;
   for (int i = 1; i <= nR; ++i) rr[i] = std::sqrt(rr[i - 1] * rr[i - 1] + r2s);

   TFile fb(baseFile);
   auto *tb = (TTree *)fb.Get("cbmsim");
   TClonesArray *eb = nullptr;
   tb->SetBranchAddress("AtEventH", &eb);
   TFile fd(dlcFile);
   auto *td = (TTree *)fd.Get("cbmsim");
   TClonesArray *ed = nullptr;
   td->SetBranchAddress("AtEventH", &ed);

   // pick a clean, well-populated baseline event
   Long64_t ie = 0;
   for (Long64_t i = 0; i < tb->GetEntries(); ++i) {
      tb->GetEntry(i);
      auto *e = (AtEvent *)(eb ? eb->At(0) : nullptr);
      if (e && e->GetNumHits() > 35 && e->GetNumHits() < 90) { ie = i; break; }
   }

   int nh; double cx, cy;
   // (1) baseline vs DLC, ZOOMED to one track -> pad-level charge spreading visible
   auto *hpB = buildPUMApads(), *hpD = buildPUMApads();
   int nhB, nhD; double bx, by, dx, dy;
   fillEvent(tb, eb, ie, hpB, nhB, bx, by);
   fillEvent(td, ed, ie, hpD, nhD, dx, dy);
   auto *c1 = new TCanvas("c1", "", 1200, 560);
   c1->Divide(2, 1);
   c1->cd(1); gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.12);
   drawZoom(hpB, nhB, bx, by, 26, "Baseline (no DLC)");
   c1->cd(2); gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.12);
   drawZoom(hpD, nhD, bx, by, 26, "DLC 1.35 M#Omega/#Box");
   c1->SaveAs(outdir + "/padplane_dlc.png");

   // (2) four single events (baseline), FULL plane: fired pads + ring circles
   auto *c2 = new TCanvas("c2", "", 1150, 1000);
   c2->Divide(2, 2);
   std::vector<TH2Poly *> hps = {buildPUMApads(), buildPUMApads(), buildPUMApads(), buildPUMApads()};
   int panel = 0;
   for (Long64_t i = 0; i < tb->GetEntries() && panel < 4; ++i) {
      tb->GetEntry(i);
      auto *e = (AtEvent *)(eb ? eb->At(0) : nullptr);
      if (!e || e->GetNumHits() < 35) continue;
      fillEvent(tb, eb, i, hps[panel], nh, cx, cy);
      c2->cd(panel + 1); gPad->SetRightMargin(0.16); gPad->SetLeftMargin(0.12);
      drawFull(hps[panel], nh, rr, Form("event %lld", i));
      panel++;
   }
   c2->SaveAs(outdir + "/padplane_event.png");
   printf("PADPLANE_DONE base event %lld, panels %d\n", ie, panel);
}
