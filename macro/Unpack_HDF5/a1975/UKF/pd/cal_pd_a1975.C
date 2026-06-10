/// @file cal_pd_a1975.C
/// @brief Energy-calibrate the 16C(p,d)15C Ex spectrum: shift so the 15C g.s. = 0.
///
/// The reconstructed g.s./0.74 doublet has the CORRECT spacing (~0.75 MeV ≈ lit
/// 0.740), so the Ex scale is right and the residual is a pure constant offset
/// (~+0.5 MeV) — a global kinematic offset (the 16C loses ~nothing in the thin H2
/// gas, so it is NOT a vertex effect). Calibration = subtract the g.s. centroid,
/// found by fitting the doublet with the separation FIXED at 0.740 MeV (stable
/// 1-position fit). Reads the cached per-track Ex from deuteron_kin.root.
///
///   root -b -q 'pd/cal_pd_a1975.C'

void cal_pd_a1975(TString cacheFile = "deuteron_kin.root")
{
   gStyle->SetOptStat(0);
   TFile *f = TFile::Open(cacheFile);
   TNtuple *t = (TNtuple *)f->Get("dk");
   float ex;
   t->SetBranchAddress("ex", &ex);
   std::vector<float> v;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (ex > -8 && ex < 18)
         v.push_back(ex);
   }
   f->Close();

   TH1F *h = new TH1F("h", "", 240, -10, 20);
   h->SetDirectory(nullptr);
   for (float e : v)
      h->Fill(e);

   // fit g.s. doublet with separation FIXED at 0.740 MeV: g.s. mean = [1]
   TF1 dg("dg", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1]-0.740)/[4])^2)", -0.7, 2.1);
   double mx = h->GetMaximum();
   dg.SetParameters(0.5 * mx, 0.5, 0.30, 0.6 * mx, 0.30);
   dg.SetParLimits(1, -0.5, 1.2);
   dg.SetParLimits(2, 0.15, 0.6);
   dg.SetParLimits(4, 0.15, 0.6);
   h->Fit(&dg, "QRN");
   double gs = dg.GetParameter(1); // g.s. centroid at the current (192 MeV) scale
   double sig_gs = dg.GetParameter(2), sig_074 = dg.GetParameter(4);
   printf("\n=== CALIBRATION ===\n");
   printf("g.s. centroid (uncalibrated) = %.3f MeV  -> constant shift = %+.3f MeV\n", gs, -gs);
   printf("doublet widths: g.s. sigma=%.3f (FWHM %.3f), 0.74 sigma=%.3f MeV\n", sig_gs, 2.3548 * sig_gs, sig_074);

   // calibrated spectrum: Ex_cal = Ex - gs
   TH1F *hc = new TH1F("hc", "^{16}C(p,d)^{15}C calibrated;E_{x}(^{15}C) [MeV];deuteron candidates", 240, -10, 20);
   hc->SetDirectory(nullptr);
   for (float e : v)
      hc->Fill(e - gs);
   hc->SetLineColor(kBlue + 1);
   hc->SetLineWidth(2);

   TCanvas *c = new TCanvas("ccal", "cal", 1050, 720);
   hc->GetXaxis()->SetRangeUser(-3, 11);
   hc->Draw("hist");
   double lev[] = {0.0, 0.740, 3.103, 4.220, 4.657, 5.866};
   const char *lbl[] = {"g.s. 1/2^{+}", "0.74", "3.10", "4.22", "4.66", "5.87"};
   double ymax = hc->GetMaximum();
   for (int i = 0; i < 6; ++i) {
      TLine *l = new TLine(lev[i], 0, lev[i], ymax);
      l->SetLineColor(kRed);
      l->SetLineStyle(2);
      l->Draw();
      TLatex *tx = new TLatex(lev[i] + 0.05, ymax * (0.97 - 0.06 * (i % 2)), lbl[i]);
      tx->SetTextColor(kRed + 1);
      tx->SetTextSize(0.026);
      tx->Draw();
   }
   c->SaveAs("pd/plots/ex_pd_calibrated.png");
   printf("saved pd/plots/ex_pd_calibrated.png  (red dashed = known 15C levels)\n");

   // report calibrated positions of the higher peaks (local maxima after shift)
   printf("\ncalibrated peak scan (smoothed local maxima > 0.5 MeV):\n");
   TH1F *hs = (TH1F *)hc->Clone("hs");
   hs->SetDirectory(nullptr);
   hs->Rebin(2);
   hs->Smooth(2);
   for (int b = 2; b < hs->GetNbinsX(); ++b) {
      double x = hs->GetBinCenter(b), y = hs->GetBinContent(b);
      if (x > 0.4 && x < 8 && y > hs->GetBinContent(b - 1) && y >= hs->GetBinContent(b + 1) && y > 0.25 * ymax)
         printf("   peak ~ %.2f MeV\n", x);
   }
}
