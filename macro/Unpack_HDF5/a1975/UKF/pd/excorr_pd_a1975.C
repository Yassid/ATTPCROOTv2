/// @file excorr_pd_a1975.C
/// @brief Flatten the Ex-vs-theta_cm dependence to sharpen the 15C(p,d) peaks.
///
/// A discrete state must have Ex independent of theta_cm, but the reconstructed
/// g.s. band is U-shaped in theta_cm (~0.65 MeV swing) — a kinematic systematic
/// that smears the integrated spectrum. This fits the g.s.-region Ex(theta_cm)
/// profile with a quadratic and subtracts it event-by-event, so every angle is
/// brought to the same Ex. Shows before/after and the g.s. doublet resolution gain.
///
///   root -b -q 'pd/excorr_pd_a1975.C'

// g.s. centroid + width from a doublet fit with the 0.74 separation FIXED (stable)
static void gsFit(TH1 *h, double &gs, double &fwhm, double &ptv)
{
   h->GetXaxis()->SetRangeUser(-2, 2.5);
   double pk = h->GetBinCenter(h->GetMaximumBin());
   h->GetXaxis()->SetRange(0, 0);
   TF1 dg("dg", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]*exp(-0.5*((x-[1]-0.740)/[4])^2)", pk - 1.1, pk + 1.1);
   double mx = h->GetMaximum();
   dg.SetParameters(0.6 * mx, pk - 0.4, 0.30, mx, 0.30);
   dg.SetParLimits(1, pk - 1.1, pk + 0.3);
   dg.SetParLimits(2, 0.10, 0.8);
   dg.SetParLimits(4, 0.10, 0.8);
   h->Fit(&dg, "QRN");
   gs = dg.GetParameter(1);
   fwhm = 2.3548 * dg.GetParameter(2);
   // data peak-to-valley of the doublet
   double p1 = h->GetBinContent(h->GetXaxis()->FindBin(gs));
   double p2 = h->GetBinContent(h->GetXaxis()->FindBin(gs + 0.74));
   double vmin = 1e9;
   for (int b = h->GetXaxis()->FindBin(gs + 0.15); b <= h->GetXaxis()->FindBin(gs + 0.6); ++b)
      vmin = std::min(vmin, h->GetBinContent(b));
   ptv = std::min(p1, p2) / std::max(1.0, vmin);
}

void excorr_pd_a1975(double exShift = -0.615, TString cacheFile = "deuteron_kin.root")
{
   gStyle->SetOptStat(0);
   TFile *f = TFile::Open(cacheFile);
   TNtuple *t = (TNtuple *)f->Get("dk");
   float ex, thcm;
   t->SetBranchAddress("ex", &ex);
   t->SetBranchAddress("thcm", &thcm);
   std::vector<float> vex, vtc;
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (ex > -8 && ex < 18 && thcm > 10 && thcm < 160) {
         vex.push_back(ex + exShift);
         vtc.push_back(thcm);
      }
   }
   f->Close();

   // Per theta_cm bin, fit the g.s. PEAK position (lower peak of the fixed-0.74
   // doublet) — clean of the 0.74 population, unlike the region mean.
   TGraph *gpk = new TGraph();
   const int NB = 13;
   double tlo = 20, thi = 150, bw = (thi - tlo) / NB;
   for (int b = 0; b < NB; ++b) {
      double tc0 = tlo + b * bw, tc1 = tc0 + bw, tcc = 0.5 * (tc0 + tc1);
      TH1F *hb = new TH1F("hbb", "", 120, -3, 5);
      hb->SetDirectory(nullptr);
      long n = 0;
      for (size_t i = 0; i < vex.size(); ++i)
         if (vtc[i] >= tc0 && vtc[i] < tc1) {
            hb->Fill(vex[i]);
            ++n;
         }
      if (n > 200) {
         double gs, fwhm, ptv;
         gsFit(hb, gs, fwhm, ptv);
         gpk->SetPoint(gpk->GetN(), tcc, gs);
      }
      delete hb;
   }
   TF1 sys("sys", "pol2", 20, 150);
   gpk->Fit(&sys, "QRN");
   double smin = 1e9;
   for (double a = 20; a <= 150; a += 1)
      smin = std::min(smin, sys.Eval(a));
   printf("g.s.-peak vs theta_cm quadratic: %.4g %+.4g*x %+.4g*x^2 ; min=%.3f\n", sys.GetParameter(0),
          sys.GetParameter(1), sys.GetParameter(2), smin);

   // before / after spectra
   TH1F *hb = new TH1F("hb", "", 200, -6, 14);
   hb->SetDirectory(nullptr);
   TH1F *ha = new TH1F("ha", "", 200, -6, 14);
   ha->SetDirectory(nullptr);
   (void)smin;
   for (size_t i = 0; i < vex.size(); ++i) {
      hb->Fill(vex[i]);
      // flatten the g.s. to 0: subtract its theta_cm-dependent position (zero-mean by calibration)
      double tc = std::min(150.0, std::max(20.0, (double)vtc[i]));
      ha->Fill(vex[i] - sys.Eval(tc));
   }
   double gb, fb, pb, ga, fa, pa;
   gsFit(hb, gb, fb, pb);
   gsFit(ha, ga, fa, pa);
   printf("\n              g.s.pos   g.s.FWHM   doublet peak/valley\n");
   printf("  BEFORE      %+6.3f    %6.3f        %.2f\n", gb, fb, pb);
   printf("  AFTER       %+6.3f    %6.3f        %.2f\n", ga, fa, pa);
   printf("  => FWHM %.0f%% better, peak/valley %.0f%% deeper\n", 100 * (fb - fa) / fb, 100 * (pa - pb) / pb);

   TCanvas *c = new TCanvas("c", "excorr", 1400, 560);
   c->Divide(2, 1);
   c->cd(1);
   gpk->SetTitle("g.s. peak position vs #theta_{cm};#theta_{cm} [deg];E_{x}^{g.s.} [MeV]");
   gpk->SetMarkerStyle(20);
   gpk->Draw("AP");
   sys.SetLineColor(kRed);
   sys.Draw("same");
   c->cd(2);
   hb->SetLineColor(kGray + 2);
   hb->SetLineWidth(2);
   hb->SetTitle("^{15}C Ex: before (grey) vs theta_{cm}-corrected (blue);E_{x}(^{15}C) [MeV];counts");
   hb->GetXaxis()->SetRangeUser(-3, 8);
   hb->Draw("hist");
   ha->SetLineColor(kBlue + 1);
   ha->SetLineWidth(2);
   ha->Draw("hist same");
   c->SaveAs("pd/plots/excorr.png");
   printf("saved pd/plots/excorr.png\n");
}
