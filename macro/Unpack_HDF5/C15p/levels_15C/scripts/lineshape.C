#pragma once
// The g.s. lineshape. A bifurcated Gaussian leaves chi2/bin ~9 in the ATTPCROOT core and its
// missing tail is what poisons every fit above Sn. Test four shapes on the g.s. peak alone,
// in a window where no resonance contributes.
//
//   G1  single Gaussian
//   G2  bifurcated Gaussian (sigma_L, sigma_R)
//   GT  bifurcated Gaussian + EXPONENTIAL TAILS, continuous and smooth at |t| = alpha:
//          f(t) = exp(-t^2/2)                 |t| <= alpha
//               = exp(alpha^2/2 - alpha|t|)   |t| >  alpha
//   GT2 same but only a right-hand tail (the side that matters, towards the resonances)
static double shapeGT(double X, double mu, double sL, double sR, double aL, double aR)
{
   double s = (X < mu) ? sL : sR, a = (X < mu) ? aL : aR;
   double t = std::fabs(X - mu) / s;
   return (t <= a) ? std::exp(-0.5 * t * t) : std::exp(0.5 * a * a - a * t);
}
double f_G1(double* x, double* p)
{ return p[0] * std::exp(-0.5 * std::pow((x[0] - p[1]) / p[2], 2)) + p[3] + p[4] * x[0]; }
double f_G2(double* x, double* p)
{ double s = (x[0] < p[1]) ? p[2] : p[3];
  return p[0] * std::exp(-0.5 * std::pow((x[0] - p[1]) / s, 2)) + p[4] + p[5] * x[0]; }
double f_GT(double* x, double* p)
{ return p[0] * shapeGT(x[0], p[1], p[2], p[3], p[4], p[5]) + p[6] + p[7] * x[0]; }
void lineshape()
{
   TFile f("dd_ex.root");
   struct C { const char* nm; double fw, lo, hi; } CC[2] =
      {{"ATTPCROOT", 1.012, -2.0, 2.2}, {"Spyral", 0.548, -1.6, 2.2}};
   for (auto& c : CC) {
      TH1D* h = (TH1D*)f.Get(c.nm); double s0 = c.fw / 2.355;
      double pk = h->GetBinContent(h->GetMaximumBin());
      printf("\n== %s   g.s. window [%.1f, %.1f]\n", c.nm, c.lo, c.hi);
      auto score = [&](TF1& g, const char* nm, int np) {
         h->Fit(&g, "RQ0NL");
         double c2 = 0; int nd = 0;
         for (int i = 1; i <= h->GetNbinsX(); i++) { double x = h->GetBinCenter(i), n = h->GetBinContent(i);
            if (x < c.lo || x > c.hi || n < 10) continue;
            double m = g.Eval(x); c2 += (n - m) * (n - m) / n; nd++; }
         printf("   %-4s  chi2/ndf %7.2f  (%d bins, %d par)", nm, c2 / std::max(1, nd - np), nd, np);
         if (std::string(nm) == "GT") printf("   alpha_L %.2f  alpha_R %.2f  FWHM %.3f/%.3f",
            g.GetParameter(4), g.GetParameter(5), 2.355 * g.GetParameter(2), 2.355 * g.GetParameter(3));
         printf("\n");
         return c2 / std::max(1, nd - np); };
      TF1 g1("g1", f_G1, c.lo, c.hi, 5);
      g1.SetParameters(pk, 0, s0, 10, 0); g1.SetParLimits(2, 0.05, 2.0); score(g1, "G1", 5);
      TF1 g2("g2", f_G2, c.lo, c.hi, 6);
      g2.SetParameters(pk, 0, s0, s0, 10, 0);
      g2.SetParLimits(2, 0.05, 2.0); g2.SetParLimits(3, 0.05, 2.0); score(g2, "G2", 6);
      TF1 gt("gt", f_GT, c.lo, c.hi, 8);
      gt.SetParameters(pk, 0, s0, s0, 1.5, 1.5, 10, 0);
      gt.SetParLimits(1, -0.8, 0.8);
      gt.SetParLimits(2, 0.05, 2.0); gt.SetParLimits(3, 0.05, 2.0);
      gt.SetParLimits(4, 0.3, 6.0);  gt.SetParLimits(5, 0.3, 6.0);
      score(gt, "GT", 8);
      // what the tail is worth far from the peak
      printf("   model value at Ex = 2.5, 4.0, 6.0 :  G2 %8.2f %8.2f %8.2f |  GT %8.2f %8.2f %8.2f\n",
             g2.Eval(2.5) - g2.GetParameter(4) - g2.GetParameter(5) * 2.5,
             g2.Eval(4.0) - g2.GetParameter(4) - g2.GetParameter(5) * 4.0,
             g2.Eval(6.0) - g2.GetParameter(4) - g2.GetParameter(5) * 6.0,
             gt.GetParameter(0) * shapeGT(2.5, gt.GetParameter(1), gt.GetParameter(2), gt.GetParameter(3), gt.GetParameter(4), gt.GetParameter(5)),
             gt.GetParameter(0) * shapeGT(4.0, gt.GetParameter(1), gt.GetParameter(2), gt.GetParameter(3), gt.GetParameter(4), gt.GetParameter(5)),
             gt.GetParameter(0) * shapeGT(6.0, gt.GetParameter(1), gt.GetParameter(2), gt.GetParameter(3), gt.GetParameter(4), gt.GetParameter(5)));
   }
}
