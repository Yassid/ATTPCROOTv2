/// @file mda_C14.C
/// @brief Multipole decomposition of the 14C(p,p') angular distributions over {L=2, L=3}.
///
/// WHY ONLY TWO MULTIPOLES. An MDA is only a measurement if the basis shapes are distinguishable
/// over the angular range actually measured. At 25-135 deg the normalised DWBA shapes correlate as
///
///        L=0    L=2    L=3    L=4
///   L=0  1.00   0.92  -0.06   0.81
///   L=2  0.92   1.00   0.30   0.96
///   L=3 -0.06   0.30   1.00   0.44
///
/// so L=0, 2 and 4 are mutually degenerate -- L=2 against L=4 at 0.96 is the same shape within the
/// errors -- while L=3 is cleanly separable from all of them. Condition numbers of the design
/// matrix: {2,3} = 7.3 and {0,2} = 7.2, both usable; {0,2,3} = 47, {2,3,4} = 39, {0,2,3,4} = 66,
/// all of which return coefficients that look quantitative and are noise. THIS CANNOT SEARCH OVER
/// ALL L. It answers one question: is there octupole strength under what is being called a 2+.
///
/// The degeneracy is set by the angular RANGE, not the statistics: L=2 and L=4 separate forward of
/// 25 deg, which is exactly where the acceptance dips to 0.48 and the analysis stops.
///
/// The basis is built at each level's OWN excitation energy, because the DWBA shape depends on the
/// Q-value of the excitation.
///
/// FIT. Non-negative least squares on two coefficients: solve the normal equations, and if either
/// coefficient comes out negative, drop it and refit the survivor alone. Errors are from a
/// parametric bootstrap over the statistical errors of the points, which also gives the
/// correlation between the two coefficients -- the number that says whether a decomposition is
/// meaningful or is two shapes trading against each other.
///
///   root -b -q 'mda_C14.C()'
namespace mda {
TGraph *rd(TString f)
{
   auto *g = new TGraph();
   std::ifstream in(f.Data());
   double a, b;
   while (in >> a >> b) if (b > 0) g->SetPoint(g->GetN(), a, b);
   return g;
}
/// Non-negative least squares for y = a2*f2 + a3*f3, weighted by 1/e^2.
struct Fit { double a2, a3, chi2; int ndf; };
Fit nnls(const std::vector<double> &y, const std::vector<double> &e,
         const std::vector<double> &f2, const std::vector<double> &f3)
{
   double s22 = 0, s33 = 0, s23 = 0, s2y = 0, s3y = 0;
   for (size_t i = 0; i < y.size(); ++i) {
      double w = 1.0 / (e[i] * e[i]);
      s22 += w * f2[i] * f2[i]; s33 += w * f3[i] * f3[i]; s23 += w * f2[i] * f3[i];
      s2y += w * f2[i] * y[i];  s3y += w * f3[i] * y[i];
   }
   double det = s22 * s33 - s23 * s23;
   double a2 = 0, a3 = 0;
   if (std::fabs(det) > 1e-30) { a2 = (s33 * s2y - s23 * s3y) / det; a3 = (s22 * s3y - s23 * s2y) / det; }
   if (a2 < 0) { a2 = 0; a3 = s33 > 0 ? s3y / s33 : 0; }
   if (a3 < 0) { a3 = 0; a2 = s22 > 0 ? s2y / s22 : 0; }
   double c2 = 0;
   for (size_t i = 0; i < y.size(); ++i)
      c2 += std::pow((y[i] - a2 * f2[i] - a3 * f3[i]) / e[i], 2);
   return {a2, a3, c2, (int)y.size() - 2};
}
} // namespace mda

void mda_C14(TString distFile = "plots/fit_angles_ps_dist_pr00.root",
             Double_t fitLo = 25, Double_t fitHi = 135, Int_t nBoot = 2000, UInt_t seed = 20260827)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TString pdir = here + "/../ptolemy/dat/";
   TFile *fd = TFile::Open(here + "/" + distFile);
   if (!fd || fd->IsZombie()) { printf("\033[1;31mno %s\033[0m\n", distFile.Data()); return; }
   if (auto *pv = (TNamed *)fd->Get("provenance")) printf("\n  data: %s\n", pv->GetTitle());

   const int NL = 5;
   const char *gn[NL] = {"lvl0", "lvl1", "lvl2", "lvl3", "lvl4"};
   const char *nm[NL] = {"6.091 1-", "6.728 3-", "7.012 2+ (blend)", "7.341 2-", "8.317 2+"};
   const char *tg[NL] = {"6091", "6728", "7012", "7341", "8317"};
   const char *ass[NL] = {"L=1 assigned", "L=3 assigned", "L=2 assigned", "unnatural parity", "L=2 assigned"};

   TRandom3 rng(seed);
   printf("\n  ================ multipole decomposition over {L=2, L=3} ================\n");
   printf("  basis at each level's own Ex, Perey potential, INELOCA1\n");
   printf("\n  %-18s %8s %8s %8s %8s %8s %8s  %s\n",
          "level", "f(L=2)", "f(L=3)", "chi2/ndf", "fall:dat", "L=2", "L=3", "SPANS?");
   printf("  %s\n", "  'fall' is the ratio of the distribution at 25 deg to 135 deg. If the data fall"
                     " much faster than\n  every basis shape, the basis cannot represent them and the"
                     " coefficients are meaningless.");

   auto *c1 = new TCanvas("cmda", "", 1500, 900); c1->Divide(3, 2);
   for (int i = 0; i < NL; ++i) {
      auto *gd = (TGraphErrors *)fd->Get(gn[i]);
      auto *g2 = mda::rd(pdir + TString::Format("mdaP_%s_L2.dat", tg[i]));
      auto *g3 = mda::rd(pdir + TString::Format("mdaP_%s_L3.dat", tg[i]));
      if (!gd || !g2->GetN() || !g3->GetN()) { printf("  %-18s missing input\n", nm[i]); continue; }

      std::vector<double> th, y, e, f2, f3;
      for (int j = 0; j < gd->GetN(); ++j) {
         double x = gd->GetX()[j], v = gd->GetY()[j], er = gd->GetEY()[j];
         if (x < fitLo || x > fitHi || v <= 0) continue;
         if (er <= 0) er = std::sqrt(std::fabs(v));            // never a zero error
         th.push_back(x); y.push_back(v); e.push_back(er);
         f2.push_back(g2->Eval(x)); f3.push_back(g3->Eval(x));
      }
      if (th.size() < 4) { printf("  %-18s too few points\n", nm[i]); continue; }

      auto best = mda::nnls(y, e, f2, f3);
      // integrated strength each component carries over the fitted range, which is what "fraction
      // of L=3" should mean -- a ratio of coefficients alone is meaningless, the shapes differ
      double i2 = 0, i3 = 0;
      for (size_t k = 0; k < th.size(); ++k) { i2 += best.a2 * f2[k]; i3 += best.a3 * f3[k]; }
      double tot = i2 + i3;

      // parametric bootstrap: resample every point within its statistical error
      std::vector<double> b2, b3;
      for (int b = 0; b < nBoot; ++b) {
         std::vector<double> yy(y.size());
         for (size_t k = 0; k < y.size(); ++k) yy[k] = rng.Gaus(y[k], e[k]);
         auto r = mda::nnls(yy, e, f2, f3);
         double j2 = 0, j3 = 0;
         for (size_t k = 0; k < th.size(); ++k) { j2 += r.a2 * f2[k]; j3 += r.a3 * f3[k]; }
         if (j2 + j3 > 0) { b2.push_back(j2 / (j2 + j3)); b3.push_back(j3 / (j2 + j3)); }
      }
      double m3 = 0, s3 = 0, m2 = 0;
      for (double v : b3) m3 += v; m3 /= std::max<size_t>(1, b3.size());
      for (double v : b2) m2 += v; m2 /= std::max<size_t>(1, b2.size());
      for (double v : b3) s3 += (v - m3) * (v - m3);
      s3 = std::sqrt(s3 / std::max<size_t>(1, b3.size()));
      // correlation between the two fractions; they sum to one, so -1 by construction is expected
      // and the useful number is the SPREAD of f(L=3), not the correlation
      double cv = 0;
      for (size_t k = 0; k < b2.size(); ++k) cv += (b2[k] - m2) * (b3[k] - m3);
      cv /= std::max<size_t>(1, b2.size());
      double sd2 = 0; for (double v : b2) sd2 += (v - m2) * (v - m2);
      sd2 = std::sqrt(sd2 / std::max<size_t>(1, b2.size()));
      double rho = (sd2 > 0 && s3 > 0) ? cv / (sd2 * s3) : 0;

      // THE SPANNING TEST. An MDA presupposes the basis can represent the data. Compare how
      // steeply each falls across the fitted range: if the data outrun every basis shape, no
      // combination with positive coefficients can reach them, the fit pins a coefficient at zero
      // and reports a decomposition that is an artefact of the basis being too flat.
      double dFall = y.front() / std::max(1e-12, y.back());
      double f2Fall = f2.front() / std::max(1e-12, f2.back());
      double f3Fall = f3.front() / std::max(1e-12, f3.back());
      bool spans = dFall <= 1.5 * std::max(f2Fall, f3Fall) && dFall >= 0.67 * std::min(f2Fall, f3Fall);
      printf("  %-18s %8.3f %5.3f+-%.3f %8.2f %8.1f %8.1f %8.1f  %s\n", nm[i],
             tot > 0 ? i2 / tot : 0, tot > 0 ? i3 / tot : 0, s3,
             best.ndf > 0 ? best.chi2 / best.ndf : 0, dFall, f2Fall, f3Fall,
             spans ? "yes" : "NO -- result meaningless");

      c1->cd(i + 1); gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
      double ym = 0, yl = 1e30;
      for (size_t k = 0; k < y.size(); ++k) { ym = std::max(ym, y[k] + e[k]); yl = std::min(yl, y[k]); }
      auto *fr = gPad->DrawFrame(20, 0.35 * yl, 142, 3.0 * ym);
      fr->SetTitle(Form("%s;#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]", nm[i]));
      auto *q2 = new TGraph(), *q3 = new TGraph(), *qs = new TGraph();
      for (int k = 0; k < g2->GetN(); ++k) {
         double x = g2->GetX()[k];
         if (x < 20 || x > 142) continue;
         double v2 = best.a2 * g2->Eval(x), v3 = best.a3 * g3->Eval(x);
         if (v2 > 0) q2->SetPoint(q2->GetN(), x, v2);
         if (v3 > 0) q3->SetPoint(q3->GetN(), x, v3);
         if (v2 + v3 > 0) qs->SetPoint(qs->GetN(), x, v2 + v3);
      }
      q2->SetLineColor(kAzure + 2); q2->SetLineWidth(2); q2->SetLineStyle(2); q2->Draw("L same");
      q3->SetLineColor(kGreen + 3); q3->SetLineWidth(2); q3->SetLineStyle(2); q3->Draw("L same");
      qs->SetLineColor(kRed + 1);   qs->SetLineWidth(3); qs->Draw("L same");
      gd->SetMarkerStyle(20); gd->SetMarkerSize(1.0); gd->SetLineWidth(2); gd->Draw("P same");
      TLatex tx; tx.SetNDC(); tx.SetTextSize(0.050);
      tx.DrawLatex(0.15, 0.26, Form("f(L=3) = %.2f #pm %.2f", tot > 0 ? i3 / tot : 0, s3));
      tx.DrawLatex(0.15, 0.19, Form("#chi^{2}/ndf = %.2f", best.ndf > 0 ? best.chi2 / best.ndf : 0));
   }
   c1->cd(6);
   TLatex t; t.SetNDC(); t.SetTextSize(0.052);
   t.DrawLatex(0.05, 0.86, "{L=2, L=3} only.");
   t.SetTextSize(0.044);
   t.DrawLatex(0.05, 0.75, "L=0, 2 and 4 are degenerate over");
   t.DrawLatex(0.05, 0.68, "25-135 deg (corr 0.81-0.96), so a");
   t.DrawLatex(0.05, 0.61, "wider basis is unstable: cond 39-66");
   t.DrawLatex(0.05, 0.54, "against 7.3 for this one.");
   t.DrawLatex(0.05, 0.42, "L=3 IS separable (corr 0.30 to L=2),");
   t.DrawLatex(0.05, 0.35, "so this tests one thing: octupole");
   t.DrawLatex(0.05, 0.28, "strength under an assigned 2+.");
   t.SetTextColor(kRed + 1); t.SetTextSize(0.040);
   t.DrawLatex(0.05, 0.14, "It cannot search over all L.");
   TString out = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";
   c1->SaveAs(out + "17_mda_L2L3.png");
   printf("\n  wrote %s17_mda_L2L3.png\n\n", out.Data());
}
