#include "globalfit.C"
// LIKELIHOOD-RATIO test on the state window [2.2, 20], which is what the fit actually minimises.
// (Comparing the models by a Pearson chi2 restricted to bins with n >= 10 is inconsistent with the
//  Poisson likelihood being optimised and made FLATEXP look worse than the FLAT model nested
//  inside it -- impossible, and purely an artefact of the statistic.)
//
// Every alternative reduces to FLAT:   LIN c1->0 | FLATEXP A->0 | PS N->0 | THR A->0
// so Delta(-2lnL) against FLAT is a chi2 with Delta(npar) degrees of freedom.
void lrtest()
{
   const char* MN[5] = {"LIN    ", "PS     ", "THR    ", "FLAT   ", "FLATEXP"};
   int NB[5] = {2, 3, 3, 1, 3};                       // background parameters
   struct { const char* nm; double kb, fw, lo; } C[2] =
      {{"ATTPCROOT", 191.25, 1.012, -1.6}, {"Spyral", 195.00, 0.548, -1.2}};
   for (auto& c : C) {
      printf("\n================ %s ================\n", c.nm);
      double nllF = 0; Out oF;
      for (int m : {3}) { oF = bestFit(c.nm, c.kb, c.fw, c.lo, m, false); nllF = oF.nll; }
      printf("bkg      bkgpar   -2lnL       vs FLAT        p (chi2, dof)   state A            state B\n");
      for (int m : {3, 0, 4, 1, 2}) {
         Out o = bestFit(c.nm, c.kb, c.fw, c.lo, m, false);
         double d = nllF - o.nll; int dof = NB[m] - 1;
         printf("%s   %d    %9.1f", MN[m], NB[m], o.nll);
         if (m == 3) printf("     --              --        ");
         else {
            double p = (d > 0 && dof > 0) ? TMath::Prob(d, dof) : 1.0;
            printf("   %+8.1f      %s", d,
                   p < 1e-4 ? "< 1e-4      " : Form("%8.4f    ", p));
         }
         printf("   %.3f(%.3f)      %.3f(%.3f)\n", o.Amu, o.dAmu, o.Bmu, o.dBmu);
      }
   }
}
