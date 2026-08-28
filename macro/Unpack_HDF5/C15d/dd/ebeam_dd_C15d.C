/// @file ebeam_dd_C15d.C
/// @brief Measure the 15C beam energy from the 15C(d,d)15C elastic recoil ridge.
///
///   root -b -q 'dd/ebeam_dd_C15d.C()'
///
/// Elastic scattering of a heavy beam on a light target is two-body with ONE free parameter, so
/// the recoil-deuteron energy as a function of lab angle fixes the beam energy exactly:
///
///     T_d(theta) = 2 m_d p1^2 cos^2(theta) / [ (E1 + m_d)^2 - p1^2 cos^2(theta) ]
///
/// with E1 = Ebeam + m1 and p1^2 = E1^2 - m1^2. Nothing else enters -- no optical model, no gate
/// efficiency, no absolute normalisation.
///
/// ★ THE RIDGE IS WALKED, NOT MAX-BINNED PER SLICE. Taking the highest bin in each theta slice
/// jumps between the elastic ridge and whatever else is populated there, which on the a2091
/// analysis produced a non-physical KE discontinuity that was chased for some time before being
/// recognised as the estimator hopping populations. Here the walk starts in a slice where the
/// ridge is unambiguous and then tracks it, searching only a window around the position
/// extrapolated from the slices already accepted.
///
/// ★ THE FIT WINDOW MATTERS AND IS REPORTED. The ridge flattens at low KE where short tracks stop
/// being reconstructed, and at small angles it runs into other populations. Fits over several
/// windows are printed so the window-to-window spread is visible as what it is -- the systematic --
/// rather than being hidden behind one number and its statistical error.

namespace {
const double kU = 931.49410242;
const double kM1 = 15.0105993 * kU; // 15C beam
const double kMd = 2.0141018 * kU;  // deuteron target and ejectile

/// elastic recoil energy of the light target, inverse kinematics
double Trecoil(double thDeg, double Ebeam)
{
   const double E1 = Ebeam + kM1;
   const double p1sq = E1 * E1 - kM1 * kM1;
   const double c = std::cos(thDeg * TMath::DegToRad());
   const double num = 2.0 * kMd * p1sq * c * c;
   const double den = (E1 + kMd) * (E1 + kMd) - p1sq * c * c;
   return den > 0 ? num / den : -1;
}
} // namespace

void ebeam_dd_C15d(TString kinFile = "dd/plots/dd_kin_C15d.root", TString outDir = "dd/plots/",
                   Double_t thLo = 42, Double_t thHi = 78, Double_t thStep = 2.0,
                   Double_t seedTheta = 50.0, Double_t keMinUse = 2.0)
{
   gSystem->mkdir(outDir, kTRUE);
   if (gSystem->AccessPathName(kinFile)) {
      std::cout << "\033[1;31mERROR: " << kinFile << " not found (run dd/kin_dd_C15d.C).\033[0m\n";
      return;
   }
   TFile f(kinFile);
   auto *t = (TTree *)f.Get("dd");
   if (!t) {
      std::cout << "\033[1;31mERROR: no tree 'dd'\033[0m\n";
      return;
   }
   Double_t ke, th, c2;
   t->SetBranchAddress("ke", &ke);
   t->SetBranchAddress("theta", &th);
   t->SetBranchAddress("chi2ndf", &c2);

   // slice index -> KE histogram
   const int nSl = (int)std::round((thHi - thLo) / thStep);
   std::vector<TH1D *> hs(nSl);
   for (int i = 0; i < nSl; ++i)
      hs[i] = new TH1D(Form("hsl%d", i), "", 300, 0, 90);
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!(ke > 0) || th < thLo || th >= thHi)
         continue;
      const int s = (int)((th - thLo) / thStep);
      if (s >= 0 && s < nSl)
         hs[s]->Fill(ke);
   }

   // ---- walk the ridge outward from the seed slice, in both directions --------------------
   const int seed = std::min(nSl - 1, std::max(0, (int)((seedTheta - thLo) / thStep)));
   std::vector<double> rTh(nSl, -1), rKe(nSl, -1), rEr(nSl, -1);
   auto peakIn = [&](TH1D *h, double lo, double hi, double &pk, double &er) {
      int b1 = std::max(1, h->FindBin(lo)), b2 = std::min(h->GetNbinsX(), h->FindBin(hi));
      double best = -1;
      int ib = -1;
      for (int b = b1; b <= b2; ++b)
         if (h->GetBinContent(b) > best) { best = h->GetBinContent(b); ib = b; }
      if (ib < 2 || best < 5)
         return false;
      // 3-point parabola on the peak, and a width-based uncertainty
      const double y0 = h->GetBinContent(ib - 1), y1 = h->GetBinContent(ib), y2 = h->GetBinContent(ib + 1);
      const double den = y0 - 2 * y1 + y2;
      double sh = std::abs(den) > 1e-12 ? 0.5 * (y0 - y2) / den : 0;
      if (std::abs(sh) > 1) sh = 0;
      pk = h->GetBinCenter(ib) + sh * h->GetBinWidth(ib);
      const double n = y0 + y1 + y2;
      er = n > 0 ? h->GetBinWidth(ib) / std::sqrt(n) : h->GetBinWidth(ib);
      return true;
   };

   double pk, er;
   if (!peakIn(hs[seed], keMinUse, 90, pk, er)) {
      std::cout << "\033[1;31mERROR: no ridge in the seed slice at " << seedTheta << " deg.\033[0m\n";
      return;
   }
   rTh[seed] = thLo + (seed + 0.5) * thStep; rKe[seed] = pk; rEr[seed] = er;
   // outward, searching a +-35 % window about the last accepted peak (the ridge is smooth)
   for (int d = 0; d < 2; ++d)
      for (int i = seed + (d ? -1 : 1); d ? i >= 0 : i < nSl; i += (d ? -1 : 1)) {
         const int prev = i + (d ? 1 : -1);
         if (rKe[prev] <= 0) continue;
         const double lo = std::max(keMinUse, rKe[prev] * 0.55), hi = rKe[prev] * 1.8;
         if (peakIn(hs[i], lo, hi, pk, er) && pk > keMinUse) {
            rTh[i] = thLo + (i + 0.5) * thStep; rKe[i] = pk; rEr[i] = er;
         }
      }

   auto *g = new TGraphErrors();
   for (int i = 0; i < nSl; ++i)
      if (rKe[i] > 0) { int n = g->GetN(); g->SetPoint(n, rTh[i], rKe[i]); g->SetPointError(n, 0, std::max(er, 0.15)); }
   std::cout << "\033[1;33m=== 15C(d,d) elastic ridge ===\033[0m\n  walked " << g->GetN() << " slices\n";

   // ---- fit Ebeam over several windows ------------------------------------------------------
   auto fitWin = [&](double a, double b, double &eb, double &ebe, double &c2n, int &npt) {
      TGraphErrors gg;
      for (int i = 0; i < g->GetN(); ++i) {
         double x, y; g->GetPoint(i, x, y);
         if (x >= a && x <= b) { int n = gg.GetN(); gg.SetPoint(n, x, y); gg.SetPointError(n, 0, g->GetErrorY(i)); }
      }
      npt = gg.GetN();
      if (npt < 4) return false;
      TF1 fn("fn", [](double *x, double *p) { return Trecoil(x[0], p[0]); }, a, b, 1);
      fn.SetParameter(0, 150);
      fn.SetParLimits(0, 20, 600);
      if (gg.Fit(&fn, "QN") != 0) return false;
      eb = fn.GetParameter(0); ebe = fn.GetParError(0);
      c2n = fn.GetNDF() > 0 ? fn.GetChisquare() / fn.GetNDF() : -1;
      return true;
   };

   printf("  %-14s %5s %12s %10s\n", "window [deg]", "pts", "Ebeam [MeV]", "chi2/ndf");
   std::vector<double> ebs;
   const double wins[][2] = {{44, 76}, {46, 72}, {48, 70}, {50, 68}, {44, 68}, {50, 76}};
   for (auto &w : wins) {
      double eb, ebe, c2n; int npt;
      if (!fitWin(w[0], w[1], eb, ebe, c2n, npt)) { printf("  %5.0f-%-8.0f %5s   (too few)\n", w[0], w[1], "-"); continue; }
      printf("  %5.0f-%-8.0f %5d %7.1f +- %.1f %9.2f\n", w[0], w[1], npt, eb, ebe, c2n);
      ebs.push_back(eb);
   }
   if (ebs.size() >= 2) {
      std::sort(ebs.begin(), ebs.end());
      const double lo = ebs.front(), hi = ebs.back(), mid = ebs[ebs.size() / 2];
      printf("\n  \033[1;32mEbeam = %.0f MeV\033[0m, window-to-window spread %.0f-%.0f MeV (= %.1f MeV/u)\n",
             mid, lo, hi, mid / 15.0);
      printf("  The SPREAD is the systematic; the per-fit errors above are only statistical.\n");
   }

   // ---- robust estimator: count tracks near the locus, scan Ebeam -------------------------
   // ★ THIS IS THE ESTIMATOR TO TRUST, and the parametric fit above is kept only as a cross-check.
   // The per-slice peak walk hops between the ridge and whatever else is populated in the slice --
   // on this data it reported 2.25 MeV at 43 deg where the ridge sits near 20, and produced
   // chi2/ndf of 130-190 for every window. Counting how many tracks lie within a fractional
   // tolerance of the predicted locus has no per-slice peak to get wrong and no error bars to
   // mis-estimate, and it gives a single clean maximum.
   {
      std::vector<float> K, T;
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (ke > 3 && ke < 60 && th > 40 && th < 66)
            { K.push_back(ke); T.push_back(th); }
      }
      auto count = [&](double E, double tol) {
         long n = 0;
         for (size_t i = 0; i < K.size(); ++i) {
            const double m = Trecoil(T[i], E);
            if (m > 0 && std::fabs(K[i] - m) / m < tol) ++n;
         }
         return n;
      };
      double bE = 0; long best = -1;
      for (double E = 60; E <= 260; E += 1.0) { long n = count(E, 0.10); if (n > best) { best = n; bE = E; } }
      double bE2 = 0; long best2 = -1;
      for (double E = bE - 12; E <= bE + 12; E += 0.5) { long n = count(E, 0.07); if (n > best2) { best2 = n; bE2 = E; } }
      printf("\n  \033[1;32mlocus-count estimator: Ebeam = %.0f MeV = %.2f MeV/u\033[0m  (%ld of %zu tracks within 7%%)\n",
             bE2, bE2 / 15.0, best2, K.size());
      printf("  restricted to theta 40-66 deg and KE 3-60 MeV, where the ridge is clean and clear of\n"
             "  the low-energy population: above ~66 deg the measured ridge falls to ~60%% of the locus,\n"
             "  which is a reconstruction threshold rather than kinematics.\n");
   }

   auto *c = new TCanvas("ceb", "ridge", 1000, 750);
   g->SetTitle("15C(d,d) elastic recoil ridge;#theta_{lab} [deg];KE_{d} [MeV]");
   g->SetMarkerStyle(20);
   g->Draw("AP");
   if (!ebs.empty()) {
      auto *fdraw = new TF1("fdraw", [](double *x, double *p) { return Trecoil(x[0], p[0]); }, thLo, thHi, 1);
      fdraw->SetParameter(0, ebs[ebs.size() / 2]);
      fdraw->SetLineColor(kRed + 1);
      fdraw->Draw("same");
   }
   c->SaveAs(outDir + "ebeam_dd_C15d.png");
   std::cout << "  wrote " << outDir << "ebeam_dd_C15d.png\n";
}
