#include "globalfit.C"
// Is the 5/2+ amplitude actually MEASURABLE, or is it degenerate with the g.s. right-hand width?
// Scan: FIX sigma_R to a grid of values, refit everything else, record the 5/2+ ratio and the
// likelihood. If the likelihood is nearly flat while the ratio swings by a large factor, the
// ratio is not determined by the data -- it is set by whatever lineshape is assumed.
static double fpk2(double* x, double* p)
{
   double X = x[0], s = (X < p[1]) ? p[2] : p[3];
   return bkg(X) + p[0] * std::exp(-0.5 * std::pow((X - p[1]) / s, 2))
        + p[4] * p[0] * std::exp(-0.5 * std::pow((X - p[1] - E52) / p[3], 2));
}
void scan(const char* nm, double kb, double fw, double lo, int mode)
{
   Out o = one(nm, kb, fw, lo, mode, false);       // sets MODE and the adopted background P0,P1,P2
   TFile f("/home/yassid/Data/C15d_statefit/dd_ex.root");
   TH1D* h = (TH1D*)f.Get(nm); h->SetDirectory(nullptr);
   printf("\n== %s   scan of the FIXED g.s. right-hand width\n", nm);
   printf("   best fit: FWHM_R %.3f  ->  5/2+ ratio %.4f +- %.4f\n", o.fwR, o.r52, o.dr52);
   printf("   FWHM_R   5/2+ ratio    5/2+ yield/g.s.   -2lnL    delta\n");
   double best = 1e30; std::vector<double> ll, rr, yy, ww;
   for (double fwR = 0.40; fwR <= 1.251; fwR += 0.05) {
      TF1 g(Form("sc_%s_%d", nm, (int)(fwR * 100)), fpk2, lo, SPLIT, 5);
      g.SetParameter(0, h->GetMaximum()); g.SetParLimits(0, 0, 1e9);
      g.SetParameter(1, o.mu);            g.SetParLimits(1, -0.8, 0.8);
      g.SetParameter(2, o.fwL / 2.355);   g.SetParLimits(2, 0.05, 1.2);
      g.FixParameter(3, fwR / 2.355);
      g.SetParameter(4, o.r52);           g.SetParLimits(4, 0, 0.9);
      TFitResultPtr r = h->Fit(&g, "RQ0NLS");
      double L = 2.0 * r->MinFcnValue();
      // 5/2+ yield relative to the g.s.: amplitudes x their own sigmas
      // integrated yield ratio: each Gaussian's amplitude x its own width. The g.s. is bifurcated,
      // so its area uses (sigma_L + sigma_R)/2 with the FITTED sigma_L, not the seed.
      double sgL = g.GetParameter(2), sgR = fwR / 2.355;
      double yr = g.GetParameter(4) * sgR / (0.5 * (sgL + sgR));
      ll.push_back(L); rr.push_back(g.GetParameter(4)); yy.push_back(yr); ww.push_back(fwR);
      best = std::min(best, L);
   }
   for (size_t i = 0; i < ll.size(); i++)
      printf("   %5.2f    %.4f        %.4f          %9.1f   %+7.1f%s\n",
             ww[i], rr[i], yy[i], ll[i], ll[i] - best, (ll[i] - best) < 1.0 ? "   <- 1 sigma" : "");
}
void degen()
{
   scan("ATTPCROOT", 191.25, 1.012, -1.6, 0);
   scan("Spyral",    195.00, 0.548, -1.2, 1);
}
