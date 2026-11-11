#include "binfit.C"
// Is there strength the (p,p') six-resonance list does not account for?
// Add ONE extra barrier BW, scan its energy from 1.5 to 11 MeV, and record the likelihood gain.
// A real missing level shows as a sharp, deep minimum in -2lnL at one energy, in BOTH codes.
static T2 TX; static double AX = 0;
double fBinX(double* x, double* p)
{
   double v = fBin(x, p);
   v += std::max(0.0, p[9]) * TX.Eval(x[0]);
   return v;
}
void extrabw()
{
   for (auto& k : CFG) {
      // reference amplitudes (free pass) for stabilisation
      REFTOT = 0; for (int i = 0; i < NB; i++) REF[i] = 0;
      { double sa[NB] = {0}, tot = 0;
        for (int ib = 0; ib < 6; ib++) { BR o = fitBin(k, ib, -2.0, 12.0, false, false);
          if (o.N < 400) continue; for (int i = 0; i < NB; i++) sa[i] += o.A[i]; tot += o.N; }
        for (int i = 0; i < NB; i++) REF[i] = sa[i]; REFTOT = tot; }

      printf("\n================ %s : scan for a 7th resonance ================\n", k.nm);
      printf("   using theta_CM bins 0+1 (the statistics), Ex [-2,12], extra BW l=1, Gamma 0.20\n");
      printf("     Ex_extra   d(-2lnL) bin0   d(-2lnL) bin1   sum\n");
      double base[2] = {0, 0};
      for (int ib = 0; ib < 2; ib++) { BR o = fitBin(k, ib, -2.0, 12.0, false, true); base[ib] = o.nll; }
      double bestE = -1, bestD = 0;
      for (double Ex = 1.6; Ex <= 11.01; Ex += 0.2) {
         double d[2] = {0, 0};
         for (int ib = 0; ib < 2; ib++) {
            // re-establish this bin's globals, then add the extra BW
            BR o = fitBin(k, ib, -2.0, 12.0, false, true);
            TX.build(Ex + GSMU, 0.20, 1, SGR, -3.0, 23.0);
            TFile f("/home/yassid/Data/C15d_statefit/dd_cmbins.root");
            TH1D* h = (TH1D*)f.Get(Form("%s_cm%d", k.nm, ib));
            h = (TH1D*)h->Clone(Form("x_%s_%d_%d", k.nm, ib, (int)(Ex * 10))); h->SetDirectory(nullptr);
            TF1 g(Form("gx_%s_%d_%d", k.nm, ib, (int)(Ex * 10)), fBinX, -2.0, 12.0, 10);
            g.SetParameter(0, o.gsA); g.SetParLimits(0, 0, 1e9);
            g.SetParameter(1, o.e5A); g.SetParLimits(1, 0, 1e9);
            g.SetParameter(2, o.cs);  g.SetParLimits(2, 0, 1e7);
            for (int i = 0; i < NB; i++) { double seed = REF[i] * o.N / REFTOT;
               g.SetParameter(3 + i, seed); g.SetParLimits(3 + i, 0, std::max(3.0 * seed, 0.4 * REF[i])); }
            double tx = TX.Eval(Ex + GSMU);
            g.SetParameter(9, tx > 0 ? 10.0 / tx : 1.0); g.SetParLimits(9, 0, 1e9);
            TFitResultPtr r = h->Fit(&g, "RQ0NLS");
            d[ib] = base[ib] - 2 * r->MinFcnValue();
         }
         double sum = d[0] + d[1];
         printf("      %5.2f      %10.1f      %10.1f   %8.1f%s\n", Ex, d[0], d[1], sum,
                sum > 60 ? "   <<<" : "");
         if (sum > bestD) { bestD = sum; bestE = Ex; }
      }
      printf("   strongest missing strength at Ex = %.2f  (d(-2lnL) = %.1f)\n", bestE, bestD);
   }
}
