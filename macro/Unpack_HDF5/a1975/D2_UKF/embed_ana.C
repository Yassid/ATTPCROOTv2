/// @file embed_ana.C -- summary of embed_dp.C: survival of correctly-reconstructed sim g.s. protons per theta_cm bin,
/// alone (mode 0/1, closure) and overlaid on real data events (mode 2 proton only, mode 3 full sim event), and the
/// mode-2 survival vs the size of the data event (AtEventH hits) it was overlaid on.
/// To compare with: data/sim-reco at vz 200-450 ~0.09 (10-20), ~0.3 (20-30), ~0.5 (30-35) (2026-09-24 window study).
void embed_ana(TString pat = "/mnt/f/a1975/embed/embed_r*.root")
{
   TChain e("e");
   e.Add(pat);
   printf("%lld rows\n", e.GetEntries());
   const int NT = 6;
   double TB[NT + 1] = {10, 15, 20, 25, 30, 35, 40};
   printf("\n thcm  |  n(sim)  mode0  mode1  mode2  mode3   (fraction found)\n");
   for (int b = 0; b < NT; ++b) {
      TString c = Form("tc0>=%g&&tc0<%g", TB[b], TB[b + 1]);
      double f[4]; int n0 = 0;
      for (int m = 0; m < 4; ++m) {
         double n = e.GetEntries(c + Form("&&mode==%d", m));
         if (m == 0) n0 = n;
         f[m] = n > 0 ? e.GetEntries(c + Form("&&mode==%d&&found", m)) / n : 0;
      }
      printf(" %2.0f-%2.0f |  %5d   %.3f  %.3f  %.3f  %.3f\n", TB[b], TB[b + 1], n0, f[0], f[1], f[2], f[3]);
   }
   double NB[7] = {0, 50, 150, 300, 600, 1200, 1e6};
   printf("\n mode 2 survival vs data-event size (all thcm 10-40) and the share of data events in each class:\n");
   double nall = e.GetEntries("mode==2");
   for (int j = 0; j < 6; ++j) {
      TString c = Form("mode==2&&nhD>=%g&&nhD<%g", NB[j], NB[j + 1]);
      double n = e.GetEntries(c);
      printf("  nhD %5.0f-%-7.0f  share %.3f   found %.3f   (n %.0f)\n", NB[j], NB[j + 1], n / nall,
             n > 0 ? e.GetEntries(c + "&&found") / n : 0, n);
   }
}
