#include "globalfit.C"
// Is the 14C+n PHASE SPACE actually required by the data, or did "floor + phase space" win only
// because it is mostly a FLAT FLOOR?  Compare five backgrounds under IDENTICAL seeding:
//   FLAT     c                            1 par   no phase space, no threshold
//   FLATEXP  c + A*exp(-k*Ex)             3 par   no phase space, no threshold
//   LIN      c0 + c1*Ex                   2 par   no phase space, no threshold
//   PS       c + N*dens(Ex)*exp(-k*Ex)    3 par   the three-body phase space, threshold at Sn
//   THR      c + A*Phi((Ex-Sn)/sig)*exp() 3 par   threshold at Sn, no phase space
// and report what fraction of the background the non-flat term actually carries.
void bkgcmp()
{
   const char* MN[5] = {"LIN    ", "PS     ", "THR    ", "FLAT   ", "FLATEXP"};
   int NPAR[5] = {2, 3, 3, 1, 3};
   struct { const char* nm; double kb, fw, lo; } C[2] =
      {{"ATTPCROOT", 191.25, 1.012, -1.6}, {"Spyral", 195.00, 0.548, -1.2}};
   for (auto& c : C) {
      printf("\n================ %s ================\n", c.nm);
      printf("bkg       npar  chi2/ndf  |  state A             state B           | bkg 2/9/18 MeV | non-flat share\n");
      double aMin = 1e9, aMax = -1e9, bMin = 1e9, bMax = -1e9;
      for (int m : {0, 1, 2, 3, 4}) {
         Out o = bestFit(c.nm, c.kb, c.fw, c.lo, m, false);
         MODE = m; KB = c.kb; SIG_RES = c.fw / 2.355;
         double b2 = bkg(2.0), b9 = bkg(9.0), b18 = bkg(18.0);
         // integrate the flat part vs the rest over Sn..20
         double sFlat = 0, sRest = 0;
         for (double X = SNUC; X <= 20.0; X += 0.1) {
            double fl = (m == 0) ? 0.0 : std::max(0.0, P0);   // LIN has no separable flat term
            sFlat += fl; sRest += std::max(0.0, bkg(X) - fl);
         }
         printf("%s   %2d   %7.2f  |  %.3f(%.3f) FW %.2f  %.3f(%.3f) FW %.2f | %4.0f %4.0f %4.0f |",
                MN[m], NPAR[m] + 11, o.chi2ndf, o.Amu, o.dAmu, o.Afw, o.Bmu, o.dBmu, o.Bfw, b2, b9, b18);
         if (m == 0 || m == 3) printf("   --\n");
         else printf("   %3.0f %%\n", 100.0 * sRest / std::max(1e-9, sFlat + sRest));
         aMin = std::min(aMin, o.Amu); aMax = std::max(aMax, o.Amu);
         bMin = std::min(bMin, o.Bmu); bMax = std::max(bMax, o.Bmu);
      }
      printf("  spread over ALL FIVE backgrounds:  state A %.3f MeV   state B %.3f MeV\n",
             aMax - aMin, bMax - bMin);
   }
}
