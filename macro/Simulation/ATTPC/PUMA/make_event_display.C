/// @file make_event_display.C
/// @brief Baseline vs DLC pad-plane event display for the SAME event, with each
///        pad/hit drawn on a CHARGE (ADC) colour scale (Z-axis palette). Left =
///        no dispersion, right = DLC 1.35 MOhm/sq. Reproduces event_display.png.
/// Run: root -b -q 'make_event_display.C(3)'    // event index 3
static void drawPanel(TTree *t, TClonesArray *evt, Long64_t ie, const char *head)
{
   t->GetEntry(ie);
   auto *ev = (AtEvent *)(evt ? evt->At(0) : nullptr);
   int nh = ev ? ev->GetNumHits() : 0;
   gPad->SetLeftMargin(0.13);
   gPad->SetBottomMargin(0.12);
   gPad->SetRightMargin(0.17);
   auto *h = new TH2F(Form("h_%s", head), Form("%s  (%d hits);x [mm];y [mm]", head, nh), 170, -128, 128, 170, -128, 128);
   h->SetZTitle("hit charge [a.u.]");
   h->GetZaxis()->SetTitleOffset(1.35);
   h->GetZaxis()->SetLabelSize(0.028);
   h->GetZaxis()->SetTitleSize(0.035);
   for (int j = 0; j < nh; j++) {
      auto &hit = ev->GetHit(j);
      auto p = hit.GetPosition();
      double q = hit.GetCharge();
      int bx = h->GetXaxis()->FindBin(p.X()), by = h->GetYaxis()->FindBin(p.Y());
      if (q > h->GetBinContent(bx, by))
         h->SetBinContent(bx, by, q);
   }
   h->SetMinimum(1e-6);
   h->Draw("COLZ");
   auto *o = new TEllipse(0, 0, 121.1, 121.1); o->SetFillStyle(0); o->SetLineColor(kGray + 2); o->Draw();
   auto *in = new TEllipse(0, 0, 62.9, 62.9); in->SetFillStyle(0); in->SetLineColor(kGray + 2); in->Draw();
   gPad->Update();
}

void make_event_display(Long64_t ie = -1, TString baseFile = "data/reco_pi_base.root",
                        TString dlcFile = "data/reco_pi_dlc.root",
                        TString out = "/Users/quantumlab/fair_install/puma_slides/figs/event_display.png")
{
   gStyle->SetOptStat(0);
   gStyle->SetTextFont(62);
   gStyle->SetLabelFont(62, "xyz");
   gStyle->SetTitleFont(62, "xyz");
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(99);

   TFile fb(baseFile);
   auto *tb = (TTree *)fb.Get("cbmsim");
   TClonesArray *eb = nullptr;
   tb->SetBranchAddress("AtEventH", &eb);
   TFile fd(dlcFile);
   auto *td = (TTree *)fd.Get("cbmsim");
   TClonesArray *ed = nullptr;
   td->SetBranchAddress("AtEventH", &ed);

   // auto-pick the first event with a clean, well-populated baseline track
   if (ie < 0) {
      for (Long64_t i = 0; i < tb->GetEntries(); i++) {
         tb->GetEntry(i);
         auto *ev = (AtEvent *)(eb ? eb->At(0) : nullptr);
         if (ev && ev->GetNumHits() > 35 && ev->GetNumHits() < 90) { ie = i; break; }
      }
      if (ie < 0) ie = 0;
   }

   auto *c = new TCanvas("c", "", 1150, 540);
   c->Divide(2, 1);
   c->cd(1); drawPanel(tb, eb, ie, "Baseline (no DLC)");
   c->cd(2); drawPanel(td, ed, ie, "DLC 1.35 M#Omega/#Box");
   c->SaveAs(out);
   printf("EVDISP_DONE event %lld\n", ie);
}
