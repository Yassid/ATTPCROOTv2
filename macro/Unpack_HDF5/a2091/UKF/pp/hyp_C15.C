/// @file hyp_C15.C
/// @brief Decide WHAT causes the Ex-vs-angle correlation of the a2091 15C(p,p') elastic peak.
///
/// The elastic Ex drifts by ~0.4 MeV across the angular range. Three hypotheses flatten it,
/// and they are NOT equivalent — they put the known 15C excited states in different places,
/// which is what discriminates them:
///
///   (A) wrong BEAM ENERGY      : scan Ebeam, KE untouched
///   (B) wrong KE SCALE         : Ebeam fixed, KE -> KE*(1-eps)  (B field / radius / sin(theta)
///                                calibration; a z-scale error acts here through Brho=BR/sin(theta))
///   (C) empirical THETA_CM correction (excorr_C15.C): flattens by construction, says nothing
///                                about the cause
///
/// For each, the macro solves for the parameter that flattens the elastic ridge and then
/// reports where the excited group lands. Known 15C levels: 6.09(1-), 6.59(0+), 6.73(3-),
/// 6.90(0-), 7.01(2+), 7.34(2-), 8.32.
///
///   root -b -q 'pp/hyp_C15.C("plots/proton_kin_5run_ukf.root",161,"ukf")'
static double omegaH(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}
static double exOfH(double m1, double m2, double m3, double m4, double Kp, double thRad, double Ke)
{
   double Et1 = Kp + m1, Et3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double a = (std::cos(thRad) * omegaH(s, m1 * m1, m2 * m2) * omegaH(u, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                 (2 * m2 * m2) +
              s + u - m2 * m2;
   return a > 0 ? std::sqrt(a) - m4 : NAN;
}

void hyp_C15(TString cache = "plots/proton_kin_5run_ukf.root", double Ebeam0 = 195.0, TString tag = "ukf",
             double chi2Cut = 5.0)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(__FILE__);
   if (!cache.BeginsWith("/"))
      cache = here + "/" + cache;
   TFile *f = TFile::Open(cache);
   if (!f || f->IsZombie()) {
      printf("\033[1;31mERROR: cannot open %s\033[0m\n", cache.Data());
      return;
   }
   TNtuple *nt = (TNtuple *)f->Get("pk");
   const double U = 931.49401, m1 = 15.0105993 * U, m2 = 1.007825 * U, m3 = 1.007825 * U, m4 = 15.0105993 * U;
   Float_t ke, th, c2;
   nt->SetBranchAddress("ke", &ke);
   nt->SetBranchAddress("theta", &th);
   nt->SetBranchAddress("chi2ndf", &c2);
   std::vector<double> vKe, vTh;
   for (Long64_t i = 0; i < nt->GetEntries(); ++i) {
      nt->GetEntry(i);
      if (c2 > chi2Cut || ke <= 0 || th < 0)
         continue;
      vKe.push_back(ke);
      vTh.push_back(th);
   }
   f->Close();

   // one pass = (g.s. mu, sigma, dEx/dtheta_lab over the elastic arc, excited-group centroid)
   auto eval = [&](double Eb, double keScale) {
      TH1D h("h", "", 360, -3, 15);
      std::vector<double> ex(vKe.size());
      for (size_t i = 0; i < vKe.size(); ++i) {
         ex[i] = exOfH(m1, m2, m3, m4, Eb, vTh[i] * TMath::DegToRad(), vKe[i] * keScale);
         if (std::isfinite(ex[i]))
            h.Fill(ex[i]);
      }
      h.GetXaxis()->SetRangeUser(-3, 3);
      double pk = h.GetBinCenter(h.GetMaximumBin());
      h.GetXaxis()->SetRange(0, 0);
      TF1 g("g", "gaus", pk - 0.7, pk + 0.7);
      h.Fit(&g, "RQ0");
      double mu = g.GetParameter(1), sg = g.GetParameter(2);
      double sx = 0, sy = 0, sxx = 0, sxy = 0;
      long n = 0;
      for (size_t i = 0; i < ex.size(); ++i) {
         if (vTh[i] < 55 || vTh[i] > 85 || !std::isfinite(ex[i]) || std::fabs(ex[i] - mu) > 3 * sg)
            continue;
         sx += vTh[i];
         sy += ex[i];
         sxx += vTh[i] * vTh[i];
         sxy += vTh[i] * ex[i];
         ++n;
      }
      double slope = n > 50 ? (n * sxy - sx * sy) / (n * sxx - sx * sx) : NAN;
      TF1 gb("gb", "gaus(0)+pol1(3)", mu + 5.0, mu + 9.5);
      // seed on the tallest bin of the group window, and keep the width physical, else the
      // fit collapses onto a single bin and reports a meaningless FWHM
      int b1 = h.GetXaxis()->FindBin(mu + 5.5), b2 = h.GetXaxis()->FindBin(mu + 9.0);
      int bTop = b1;
      for (int b = b1; b <= b2; ++b)
         if (h.GetBinContent(b) > h.GetBinContent(bTop))
            bTop = b;
      gb.SetParameters(h.GetBinContent(bTop), h.GetBinCenter(bTop), 0.45, 2, 0);
      gb.SetParLimits(0, 0.2 * h.GetBinContent(bTop), 5 * h.GetBinContent(bTop));
      gb.SetParLimits(1, mu + 5.5, mu + 9.0);
      gb.SetParLimits(2, 0.15, 1.2);
      h.Fit(&gb, "RQ0");
      struct R {
         double mu, sg, slope, bump, bsg;
      };
      return R{mu, sg, slope, gb.GetParameter(1) - mu, std::fabs(gb.GetParameter(2))};
   };

   auto solve = [&](bool beamMode) { // secant on the parameter that zeroes the slope
      double x1 = beamMode ? Ebeam0 : 1.0, x2 = beamMode ? Ebeam0 * 1.1 : 0.97;
      double s1 = beamMode ? eval(x1, 1.0).slope : eval(Ebeam0, x1).slope;
      double s2 = beamMode ? eval(x2, 1.0).slope : eval(Ebeam0, x2).slope;
      for (int it = 0; it < 25; ++it) {
         if (!std::isfinite(s1) || !std::isfinite(s2) || std::fabs(s2 - s1) < 1e-12)
            break;
         double x3 = x2 - s2 * (x2 - x1) / (s2 - s1);
         if (!std::isfinite(x3))
            break;
         x1 = x2;
         s1 = s2;
         x2 = x3;
         s2 = beamMode ? eval(x2, 1.0).slope : eval(Ebeam0, x2).slope;
         if (std::fabs(s2) < 1e-6)
            break;
      }
      return x2;
   };

   printf("\n\033[1;33m===== %s : what flattens Ex(theta)? (Ebeam0 = %.1f MeV) =====\033[0m\n", tag.Data(), Ebeam0);
   auto r0 = eval(Ebeam0, 1.0);
   printf("(0) as measured        : mu %+.3f  FWHM %.3f  slope %+.5f  group at Ex=%.3f (FWHM %.2f)\n", r0.mu,
          2.355 * r0.sg, r0.slope, r0.bump, 2.355 * r0.bsg);

   double ebSol = solve(true);
   auto rA = eval(ebSol, 1.0);
   printf("(A) beam energy   %7.2f MeV (%.3f MeV/u) : mu %+.3f  FWHM %.3f  slope %+.6f  "
          "\033[1;36mgroup at Ex=%.3f\033[0m (FWHM %.2f)\n",
          ebSol, ebSol / 14.0, rA.mu, 2.355 * rA.sg, rA.slope, rA.bump, 2.355 * rA.bsg);

   double ksSol = solve(false);
   auto rB = eval(Ebeam0, ksSol);
   printf("(B) KE scale      %7.4f (%+.2f%%)          : mu %+.3f  FWHM %.3f  slope %+.6f  "
          "\033[1;36mgroup at Ex=%.3f\033[0m (FWHM %.2f)\n",
          ksSol, 100 * (ksSol - 1), rB.mu, 2.355 * rB.sg, rB.slope, rB.bump, 2.355 * rB.bsg);

   printf("    known 15C levels near the group: 6.09(1-) 6.59(0+) 6.73(3-) 6.90(0-) 7.01(2+) 7.34(2-)\n");
}
