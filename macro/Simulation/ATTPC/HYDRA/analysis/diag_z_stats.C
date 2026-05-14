/// @file diag_z_stats.C
/// @brief Histogram pad_z and UKF_z residuals relative to MC first-point z
/// over the full event sample at p=800 MeV/c.

void diag_z_stats(int P_MeV = 800, const char *runDir = "data",
                  const char *outPng = "data/z_residual.png")
{
   TString simF = Form("%s/HYDRAsim_p%d.root", runDir, P_MeV);
   TString digiF = Form("%s/output_digi_p%d.root", runDir, P_MeV);
   TString ukfF  = Form("%s/output_ukf_HYDRA_p%d.root", runDir, P_MeV);

   TFile fS(simF), fD(digiF), fU(ukfF);
   auto *tS = (TTree *)fS.Get("cbmsim");
   auto *tD = (TTree *)fD.Get("cbmsim");
   auto *tU = (TTree *)fU.Get("cbmsim");
   auto *trks = new TClonesArray("AtMCTrack");
   auto *pts  = new TClonesArray("AtMCPoint");
   auto *ev   = new TClonesArray("AtEvent");
   auto *te   = new TClonesArray("AtTrackingEvent");
   tS->SetBranchAddress("MCTrack", &trks);
   tS->SetBranchAddress("AtTpcPoint", &pts);
   tD->SetBranchAddress("AtEventH", &ev);
   tU->SetBranchAddress("AtTrackingEvent", &te);

   auto *hPad = new TH1F("hPad", ";z_{pad}^{nearest MC-first}-z_{MC}^{first} (mm);events", 80, -40, 40);
   auto *hUkf = new TH1F("hUkf", ";z_{UKF}^{POCA}-z_{MC}^{vertex}      (mm);events", 80, -40, 40);
   auto *hUkfSm = new TH1F("hUkfSm", ";z_{UKF}^{nearest MC-first}-z_{MC}^{first} (mm);events", 80, -40, 40);

   Long64_t n = std::min({tS->GetEntries(), tD->GetEntries(), tU->GetEntries()});
   int nFilled = 0;
   for (Long64_t i = 0; i < n; ++i) {
      tS->GetEntry(i); tD->GetEntry(i); tU->GetEntry(i);
      if (trks->GetEntries() == 0) continue;
      auto *mc = (AtMCTrack *)trks->At(0);
      double mcVz = mc->GetStartZ() * 10.;

      double xFirst = 1e9, yFirst = 0, zFirst = 0;
      for (int j = 0; j < pts->GetEntries(); ++j) {
         auto *p = (AtMCPoint *)pts->At(j);
         if (p->GetVolName() != TString("drift_volume")) continue;
         if (p->GetTrackID() != 0) continue;
         double xm = p->GetX()*10, ym = p->GetY()*10, zm = p->GetZ()*10;
         if (xm < xFirst) { xFirst = xm; yFirst = ym; zFirst = zm; }
      }
      if (xFirst > 1e8) continue;

      if (ev->GetEntries() > 0) {
         double dMin = 1e9, zPad = 0;
         for (auto &h : ((AtEvent *)ev->At(0))->GetHits()) {
            auto pos = h->GetPosition();
            double d2 = (pos.X()-xFirst)*(pos.X()-xFirst) + (pos.Y()-yFirst)*(pos.Y()-yFirst);
            if (d2 < dMin) { dMin = d2; zPad = pos.Z(); }
         }
         if (dMin < 4) hPad->Fill(zPad - zFirst); // require <2mm xy distance
      }

      if (te->GetEntries() > 0) {
         auto *trkEvt = (AtTrackingEvent *)te->At(0);
         auto &fitted = trkEvt->GetFittedTracks();
         AtFittedTrack *best = nullptr; double bestChi = 1e30;
         for (auto &t : fitted) {
            if (!t->GetTrackMetadata()) continue;
            double ndf = t->GetTrackMetadata()->GetNdf();
            double cc = ndf > 0 ? t->GetTrackMetadata()->GetChi2() / ndf : 1e30;
            if (cc < bestChi) { bestChi = cc; best = t.get(); }
         }
         if (best) {
            auto v = best->GetVertex();
            hUkf->Fill(v.Z() - mcVz);
            double dMin = 1e9, zU = 0;
            for (auto &pp : best->GetSmoothedPositions()) {
               double d2 = (pp.X()-xFirst)*(pp.X()-xFirst) + (pp.Y()-yFirst)*(pp.Y()-yFirst);
               if (d2 < dMin) { dMin = d2; zU = pp.Z(); }
            }
            if (dMin < 4) hUkfSm->Fill(zU - zFirst);
            nFilled++;
         }
      }
   }

   auto *c = new TCanvas("c", "z residuals", 1500, 500);
   c->Divide(3, 1, 0.012, 0.012);
   auto fitGauss = [](TH1F *h, const char *name) {
      h->Fit("gaus", "Q0");
      auto *f = h->GetFunction("gaus");
      h->Draw();
      if (f) f->Draw("same");
      auto *l = new TLatex();
      l->SetNDC(true); l->SetTextSize(0.05);
      if (f) l->DrawLatex(0.55, 0.85, Form("#mu=%.2f mm", f->GetParameter(1)));
      if (f) l->DrawLatex(0.55, 0.78, Form("#sigma=%.2f mm", f->GetParameter(2)));
      l->DrawLatex(0.55, 0.71, Form("N=%d", (int)h->GetEntries()));
   };
   c->cd(1); fitGauss(hPad, "Pad");
   c->cd(2); fitGauss(hUkf, "UKF POCA");
   c->cd(3); fitGauss(hUkfSm, "UKF smoothed");
   c->SaveAs(outPng);
   std::cout << "Wrote " << outPng << " (N=" << nFilled << " UKF fits)\n";
}
