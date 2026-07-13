/// @file display_tracks_charge.C
/// @brief Event displays of reconstructed PUMA tracks with each pad/hit drawn
///        on a CHARGE colour scale (Z-axis palette = deposited charge/ADC).
/// Run: root -b -q 'display_tracks_charge.C("data/reco_pid_base.root")'
void display_tracks_charge(TString file = "data/reco_pid_base.root",
                           TString out = "/Users/quantumlab/fair_install/puma_slides/figs/track_charge.png")
{
   gStyle->SetOptStat(0);
   gStyle->SetTextFont(62);
   gStyle->SetLabelFont(62, "xyz");
   gStyle->SetTitleFont(62, "xyz");
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(99);

   TFile f(file);
   auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *evt = nullptr;
   t->SetBranchAddress("AtEventH", &evt);

   // pick the first 4 events that have a decent number of hits
   std::vector<Long64_t> good;
   for (Long64_t i = 0; i < t->GetEntries() && good.size() < 4; i++) {
      t->GetEntry(i);
      auto *ev = (AtEvent *)(evt ? evt->At(0) : nullptr);
      if (ev && ev->GetNumHits() > 30)
         good.push_back(i);
   }

   auto *c = new TCanvas("c", "", 1150, 1000);
   c->Divide(2, 2);
   for (size_t k = 0; k < good.size(); k++) {
      c->cd(k + 1);
      gPad->SetLeftMargin(0.13);
      gPad->SetBottomMargin(0.12);
      gPad->SetRightMargin(0.18); // room for the Z colour bar
      t->GetEntry(good[k]);
      auto *ev = (AtEvent *)evt->At(0);
      int nh = ev->GetNumHits();

      // charge map: fine bins over the pad plane, each bin = max hit charge there
      auto *h = new TH2F(Form("h%zu", k),
                         Form("event %lld  (%d hits);x [mm];y [mm]", good[k], nh),
                         170, -128, 128, 170, -128, 128);
      h->SetZTitle("hit charge [a.u.]");
      h->GetZaxis()->SetTitleOffset(1.35);
      h->GetZaxis()->SetLabelSize(0.028);
      h->GetZaxis()->SetTitleSize(0.035);
      for (int j = 0; j < nh; j++) {
         auto &hit = ev->GetHit(j);
         auto p = hit.GetPosition();
         double q = hit.GetCharge();
         int bx = h->GetXaxis()->FindBin(p.X());
         int by = h->GetYaxis()->FindBin(p.Y());
         if (q > h->GetBinContent(bx, by))
            h->SetBinContent(bx, by, q);
      }
      h->SetMinimum(1e-6); // leave empty bins white
      h->Draw("COLZ");

      // annular PUMA pad plane outline
      auto *o = new TEllipse(0, 0, 121.1, 121.1);
      o->SetFillStyle(0);
      o->SetLineColor(kGray + 2);
      o->Draw();
      auto *in = new TEllipse(0, 0, 62.9, 62.9);
      in->SetFillStyle(0);
      in->SetLineColor(kGray + 2);
      in->Draw();
      gPad->Update();
   }
   c->SaveAs(out);
   printf("TRACKS_QMAP_DONE events %zu\n", good.size());
}
