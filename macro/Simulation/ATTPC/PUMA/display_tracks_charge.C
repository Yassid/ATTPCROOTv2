/// @file display_tracks_charge.C
/// @brief Event displays of reconstructed PUMA tracks, colour-coded by charge
///        sign (from curvature sense in the 4 T field): q>0 red, q<0 blue.
/// Run: root -b -q 'display_tracks_charge.C("data/reco_pid_base.root")'
int chargeSignOf(AtTrack &tr)
{
   const auto &hits = tr.GetHitArray();
   if (hits.size() < 3)
      return 0;
   std::vector<std::pair<double, std::pair<double, double>>> byR;
   for (const auto &h : hits) {
      const auto &p = h->GetPosition();
      byR.emplace_back(p.X() * p.X() + p.Y() * p.Y(), std::make_pair(p.X(), p.Y()));
   }
   std::sort(byR.begin(), byR.end());
   const auto &p0 = byR.front().second;
   const auto &pm = byR[byR.size() / 2].second;
   const auto &pN = byR.back().second;
   double crossZ = (pm.first - p0.first) * (pN.second - p0.second) - (pm.second - p0.second) * (pN.first - p0.first);
   return (crossZ > 0) ? -1 : +1; // calibrated: K+ (q>0) -> crossZ<0
}

void display_tracks_charge(TString file = "data/reco_pid_base.root",
                           TString out = "/Users/quantumlab/fair_install/puma_slides/figs/track_charge.png")
{
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62, "xyz"); gStyle->SetTitleFont(62, "xyz");
   gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
   TFile f(file); auto *t = (TTree *)f.Get("cbmsim"); TClonesArray *pat = nullptr;
   t->SetBranchAddress("AtPatternEvent", &pat);
   std::vector<Long64_t> good;
   for (Long64_t i = 0; i < t->GetEntries() && good.size() < 4; i++) {
      t->GetEntry(i);
      auto *pe = (AtPatternEvent *)(pat ? pat->At(0) : nullptr);
      if (pe && pe->GetTrackCand().size() >= 2) good.push_back(i);
   }
   auto *c = new TCanvas("c", "", 1100, 1000); c->Divide(2, 2);
   for (size_t k = 0; k < good.size(); k++) {
      c->cd(k + 1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.12);
      t->GetEntry(good[k]);
      auto *pe = (AtPatternEvent *)pat->At(0);
      auto *fr = new TH2F(Form("fr%zu", k), Form("event %lld  (%zu tracks);x [mm];y [mm]", good[k], pe->GetTrackCand().size()),
                          10, -130, 130, 10, -130, 130);
      fr->Draw();
      auto *o = new TEllipse(0, 0, 121.1, 121.1); o->SetFillStyle(0); o->SetLineColor(kGray + 1); o->Draw();
      auto *in = new TEllipse(0, 0, 62.9, 62.9); in->SetFillStyle(0); in->SetLineColor(kGray + 1); in->Draw();
      for (auto &tr : pe->GetTrackCand()) {
         int q = chargeSignOf(tr);
         Color_t col = (q > 0) ? kRed + 1 : (q < 0 ? kAzure + 2 : kGray + 1);
         auto *g = new TGraph();
         for (const auto &h : tr.GetHitArray()) { const auto &p = h->GetPosition(); g->SetPoint(g->GetN(), p.X(), p.Y()); }
         g->SetMarkerStyle(20); g->SetMarkerSize(0.5); g->SetMarkerColor(col); g->Draw("P same");
      }
   }
   c->cd(1);
   auto *lg = new TLegend(0.62, 0.78, 0.88, 0.9); lg->SetTextFont(62);
   auto *rp = new TMarker(0, 0, 20); rp->SetMarkerColor(kRed + 1);
   auto *bp = new TMarker(0, 0, 20); bp->SetMarkerColor(kAzure + 2);
   lg->AddEntry(rp, "q > 0", "p"); lg->AddEntry(bp, "q < 0", "p"); lg->Draw();
   c->SaveAs(out);
   printf("TRACKS_DONE events %zu\n", good.size());
}
