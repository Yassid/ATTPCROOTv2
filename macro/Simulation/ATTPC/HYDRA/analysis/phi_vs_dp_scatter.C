/// @file phi_vs_dp_scatter.C
/// @brief 2x3 grid of Δφ-vs-Δp/p scatter plots, one per momentum scan point.
/// Tests the hypothesis that σ_φ at the first cluster is driven by σ_p/p
/// (curvature uncertainty) at low p, and saturates at a pad-resolution
/// spatial floor at high p.
///
/// MC reference uses the first AtMCPoint inside drift_volume — same as
/// make_performance.C — so the apples-to-apples comparison is preserved.
///
/// Output: data/phi_vs_dp_HYDRA.png

void phi_vs_dp_scatter(const char *outPng = "data/phi_vs_dp_HYDRA.png")
{
   gStyle->SetOptStat(0);
   gStyle->SetTitleSize(0.06, "T");
   gStyle->SetTitleSize(0.055, "XY");
   gStyle->SetLabelSize(0.045, "XY");
   gStyle->SetPadLeftMargin(0.145);
   gStyle->SetPadBottomMargin(0.135);
   gStyle->SetPadRightMargin(0.04);
   gStyle->SetPadTopMargin(0.085);

   const std::vector<int> plist = {200, 400, 600, 800, 1000, 1200};
   const int nP = plist.size();
   const double mass_pi = 139.57039;
   const double radToMrad = 1000.;

   auto *c = new TCanvas("c", "Δφ vs Δp/p", 1800, 1000);
   c->Divide(3, 2, 0.012, 0.015);

   for (int k = 0; k < nP; ++k) {
      int P = plist[k];
      TFile fS(Form("data/HYDRAsim_p%d.root", P));
      TFile fU(Form("data/output_ukf_HYDRA_p%d.root", P));
      auto *tS = (TTree *)fS.Get("cbmsim");
      auto *tU = (TTree *)fU.Get("cbmsim");
      if (!tS || !tU) continue;

      auto *trks = new TClonesArray("AtMCTrack");
      auto *pts  = new TClonesArray("AtMCPoint");
      auto *te = new TClonesArray("AtTrackingEvent");
      tS->SetBranchAddress("MCTrack", &trks);
      tS->SetBranchAddress("AtTpcPoint", &pts);
      tU->SetBranchAddress("AtTrackingEvent", &te);

      // 2D histogram — wider Δp/p axis for low-p (more spread)
      double dpRange = (P <= 300) ? 0.30 : 0.20;
      double dphiRange = (P <= 300) ? 250. : 100.;
      auto *h2 = new TH2F(Form("h2_%d", P),
                          Form("p = %d MeV/c;#Deltap / p_{MC};#Delta#varphi (mrad)", P),
                          80, -dpRange, dpRange, 80, -dphiRange, dphiRange);
      h2->SetDirectory(nullptr);

      Long64_t n = std::min(tS->GetEntries(), tU->GetEntries());
      for (Long64_t i = 0; i < n; ++i) {
         tS->GetEntry(i);
         tU->GetEntry(i);
         if (trks->GetEntries() == 0) continue;
         auto *mc = (AtMCTrack *)trks->At(0);
         if (std::abs(mc->GetPdgCode()) != 211) continue;
         double Px = mc->GetPx(), Py = mc->GetPy(), Pz = mc->GetPz();
         double pmcGeV = std::sqrt(Px * Px + Py * Py + Pz * Pz);
         double pmc = pmcGeV * 1000.;

         // Angles at first MC point in drift_volume
         double phMC = std::atan2(Py, Px);
         double xMinPt = 1e9;
         for (int j = 0; j < pts->GetEntries(); ++j) {
            auto *p = (AtMCPoint *)pts->At(j);
            if (p->GetVolName() != TString("drift_volume")) continue;
            if (p->GetTrackID() != 0) continue;
            double x_mm = p->GetX() * 10.;
            if (x_mm < xMinPt) {
               xMinPt = x_mm;
               double pyp = p->GetPy(), pxp = p->GetPx();
               phMC = std::atan2(pyp, pxp);
            }
         }

         if (te->GetEntries() == 0) continue;
         auto *trkEvt = (AtTrackingEvent *)te->At(0);
         auto &fitted = trkEvt->GetFittedTracks();
         AtFittedTrack *best = nullptr;
         double bestChi = 1e30;
         for (auto &t : fitted) {
            if (!t->GetTrackMetadata()) continue;
            double ndf = t->GetTrackMetadata()->GetNdf();
            double cc = ndf > 0 ? t->GetTrackMetadata()->GetChi2() / ndf : 1e30;
            if (cc < bestChi) { bestChi = cc; best = t.get(); }
         }
         if (!best) continue;
         auto kin = best->GetKinematics();
         double Efit = kin.kineticEnergy + mass_pi;
         double pfit = std::sqrt(Efit * Efit - mass_pi * mass_pi);
         double dp = (pfit - pmc) / pmc;
         double dPhi = kin.phi - phMC;
         while (dPhi > M_PI) dPhi -= 2 * M_PI;
         while (dPhi < -M_PI) dPhi += 2 * M_PI;
         h2->Fill(dp, dPhi * radToMrad);
      }

      // Linear fit y = a + b*x
      double sumW = 0, sumX = 0, sumY = 0, sumXX = 0, sumXY = 0, sumYY = 0;
      for (int ix = 1; ix <= h2->GetNbinsX(); ++ix) {
         for (int iy = 1; iy <= h2->GetNbinsY(); ++iy) {
            double w = h2->GetBinContent(ix, iy);
            if (w <= 0) continue;
            double xc = h2->GetXaxis()->GetBinCenter(ix);
            double yc = h2->GetYaxis()->GetBinCenter(iy);
            sumW += w; sumX += w * xc; sumY += w * yc;
            sumXX += w * xc * xc; sumXY += w * xc * yc; sumYY += w * yc * yc;
         }
      }
      double xbar = sumX / sumW, ybar = sumY / sumW;
      double Sxx = sumXX - sumW * xbar * xbar;
      double Sxy = sumXY - sumW * xbar * ybar;
      double Syy = sumYY - sumW * ybar * ybar;
      double slope = Sxx > 0 ? Sxy / Sxx : 0;
      double intercept = ybar - slope * xbar;
      double rho = (Sxx > 0 && Syy > 0) ? Sxy / std::sqrt(Sxx * Syy) : 0;

      c->cd(k + 1);
      gPad->SetGrid();
      gPad->SetLogz();
      h2->Draw("colz");

      // Linear-fit overlay
      auto *fline = new TF1(Form("fline_%d", P), "[0] + [1]*x", -dpRange, dpRange);
      fline->SetParameters(intercept, slope);
      fline->SetLineColor(kRed + 1);
      fline->SetLineWidth(2);
      fline->Draw("same");

      // Annotation
      auto *txt = new TLatex();
      txt->SetTextSize(0.05);
      txt->SetTextColor(kRed + 1);
      txt->SetNDC();
      txt->DrawLatex(0.18, 0.88, Form("slope = %.0f mrad", slope));
      txt->DrawLatex(0.18, 0.82, Form("#rho = %.2f", rho));
   }

   c->SaveAs(outPng);
   std::cout << "Wrote " << outPng << "\n";
}
