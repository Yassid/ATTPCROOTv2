/// @file fit_examples.C
/// @brief Gallery of Kalman-fit examples: per event, the reconstructed hits (grey)
///        with the UKF (blue) and GENFIT (red) fitted trajectories (smoothed states)
///        and their back-extrapolated vertices (stars). Shows the fit quality.
/// Run: root -b -q 'fit_examples.C("data/hs_reco375.root",12)'
void fit_examples(TString file = "data/hs_reco375.root", int nShow = 12,
                  TString out = "/Users/quantumlab/fair_install/puma_slides/figs/fit_examples.png")
{
   gSystem->Load("libAtReconstruction.so");
   gStyle->SetOptStat(0); gStyle->SetTextFont(62); gStyle->SetLabelFont(62, "xyz"); gStyle->SetTitleFont(62, "xyz");
   TFile f(file); auto *t = (TTree *)f.Get("cbmsim");
   TClonesArray *pat = nullptr, *ukf = nullptr, *gf = nullptr;
   t->SetBranchAddress("AtPatternEvent", &pat);
   t->SetBranchAddress("AtTrackingEventUKF", &ukf);
   t->SetBranchAddress("AtTrackingEventGenfit", &gf);

   std::vector<Long64_t> good;
   for (Long64_t i = 0; i < t->GetEntries() && (int)good.size() < nShow; ++i) { t->GetEntry(i);
      auto *pe = (AtPatternEvent *)(pat ? pat->At(0) : nullptr); if (!pe) continue;
      if (pe->GetTrackCand().size() < 2 || pe->GetTrackCand().size() > 4) continue;
      auto *tu = ukf->GetEntries() ? (AtTrackingEvent *)ukf->At(0) : nullptr;
      if (!tu || tu->GetFittedTracks().size() < 2) continue;
      good.push_back(i); }
   int ng = good.size(), nc = 4, nr = (ng + nc - 1) / nc;
   auto *c = new TCanvas("cfit", "", 330 * nc, 330 * nr); c->Divide(nc, nr, 0.001, 0.012);

   auto drawFits = [](TClonesArray *arr, Color_t col) {
      if (arr->GetEntries() == 0) return;
      for (const auto &ft : ((AtTrackingEvent *)arr->At(0))->GetFittedTracks()) {
         const auto &sp = ft->GetSmoothedPositions();
         if (sp.size() >= 2) { auto *pl = new TPolyLine((int)sp.size());
            for (size_t i = 0; i < sp.size(); ++i) pl->SetPoint((int)i, sp[i].X(), sp[i].Y());
            pl->SetLineColor(col); pl->SetLineWidth(2); pl->Draw(); }
         const auto &pr = ft->GetTrackPropertiesStruct();
         double vx = pr.initialPositionXtr.X(), vy = pr.initialPositionXtr.Y();
         if (std::abs(vx) > 1e-9 || std::abs(vy) > 1e-9) { auto *m = new TMarker(vx, vy, 29);
            m->SetMarkerColor(col); m->SetMarkerSize(1.8); m->Draw(); } } };

   for (int k = 0; k < ng; ++k) {
      c->cd(k + 1); gPad->SetMargin(0.03, 0.03, 0.03, 0.03); gPad->SetFixedAspectRatio(kTRUE);
      t->GetEntry(good[k]); auto *pe = (AtPatternEvent *)pat->At(0);
      auto *fr = new TH2F(Form("f%d", k), "", 1, -136, 136, 1, -136, 136);
      fr->GetXaxis()->SetLabelSize(0); fr->GetYaxis()->SetLabelSize(0);
      fr->GetXaxis()->SetTickLength(0); fr->GetYaxis()->SetTickLength(0); fr->Draw();
      for (double R : {62.9, 121.1}) { auto *e = new TEllipse(0, 0, R, R); e->SetFillStyle(0); e->SetLineColor(17); e->Draw(); }
      // reconstructed hits (grey)
      auto *gh = new TGraph();
      for (auto &tr : pe->GetTrackCand()) for (auto &h : tr.GetHitArray()) { const auto &p = h->GetPosition(); gh->SetPoint(gh->GetN(), p.X(), p.Y()); }
      gh->SetMarkerStyle(20); gh->SetMarkerSize(0.4); gh->SetMarkerColor(kGray + 1); gh->Draw("P same");
      // fitted trajectories
      drawFits(gf, kRed + 1);       // GENFIT
      drawFits(ukf, kAzure + 2);    // UKF
      auto *lab = new TLatex(0, 127, Form("event %lld", good[k]));
      lab->SetTextFont(62); lab->SetTextSize(0.056); lab->SetTextAlign(22); lab->Draw();
   }
   // legend in the first panel
   c->cd(1);
   auto *lg = new TLegend(0.02, 0.02, 0.5, 0.16); lg->SetTextFont(62); lg->SetTextSize(0.05); lg->SetBorderSize(0); lg->SetFillStyle(0);
   auto *lu = new TMarker(0, 0, 20); lu->SetMarkerColor(kAzure + 2); auto *lgf = new TMarker(0, 0, 20); lgf->SetMarkerColor(kRed + 1);
   lg->AddEntry(lu, "UKF fit", "l"); lg->AddEntry(lgf, "GENFIT fit", "l"); lg->Draw();
   c->SaveAs(out);
   printf("FIT_EXAMPLES_DONE %d events\n", ng);
}
