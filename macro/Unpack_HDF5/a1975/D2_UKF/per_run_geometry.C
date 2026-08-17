/// @file per_run_geometry.C
/// @brief Is the per-run geometry (dv / z-scale) actually varying, or is that hypothesis dead?
///
/// THE ARGUMENT THIS TESTS. Ex is computed from (KE, theta) against two-body kinematics, so it
/// IS the residual of the data from the kinematic locus. A wrong drift velocity mis-scales z,
/// hence theta, hence Ex. So the per-run Ex scale is a direct per-run geometry probe -- and it
/// does not need the CATIMA production at all.
///
/// Crucially it is measured on the matFX-OFF cache, which collapses only 0.3% of its fits in
/// every run. That sample is unbiased by the very effect under investigation, so it can say
/// whether the runs that collapse under material effects have anything geometrically wrong with
/// them in the first place.
///
///   run flat in Ex  -> geometry is uniform, the dv hypothesis is DEAD and the collapse is
///                      caused by something else entirely.
///   Ex tracks the collapse rate -> geometry really does vary and dv is worth scanning.
///
///   root -b -q per_run_geometry.C

#include <cstdio>
#include <vector>

void per_run_geometry()
{
   gStyle->SetOptStat(0);
   auto *fO = TFile::Open("/mnt/f/a1975/caches/dt_kin_dv1104.root");
   auto *fN = TFile::Open("/mnt/f/a1975/caches/dt_kin_catima.root");
   auto *tO = (TTree *)fO->Get("pk");
   auto *tN = (TTree *)fN->Get("pk");
   const char *conv = "chi2ndf<1e9";

   int runs[] = {16, 17, 18, 19, 20, 21, 22, 23, 26, 27, 31, 32, 34, 36, 37, 38, 39, 40, 41, 42,
                 43, 44, 46, 48, 57, 58, 76, 77, 78, 79, 80, 82, 83, 84, 85, 86, 87, 88, 89, 91,
                 92, 95, 96, 97, 98, 102, 103};

   std::vector<double> col, exPk, exMed, vz;
   printf("\n%6s %9s %12s %12s %12s\n", "run", "collapse", "Ex peak", "Ex median", "vertex z");
   for (int r : runs) {
      const long n = tN->GetEntries(Form("run==%d", r));
      if (n < 50) continue;
      const long b = tN->GetEntries(Form("run==%d && chi2ndf>=1e9", r));

      // matFX-OFF sample: unbiased by the collapse under study
      auto *h = new TH1D("h", "", 60, -3, 3);
      tO->Draw("ex>>h", Form("run==%d && %s", r, conv), "goff");
      if (h->GetEntries() < 100) { delete h; continue; }
      const double peak = h->GetBinCenter(h->GetMaximumBin());
      double q = 0.5, med = 0;
      h->GetQuantiles(1, &med, &q);

      auto *hz = new TH1D("hz", "", 100, 0, 1000);
      tN->Draw("vertexz>>hz", Form("run==%d && %s", r, conv), "goff");

      col.push_back(100.0 * b / n);
      exPk.push_back(peak);
      exMed.push_back(med);
      vz.push_back(hz->GetMean());
      printf("%6d %8.1f%% %12.3f %12.3f %12.1f\n", r, 100.0 * b / n, peak, med, hz->GetMean());
      delete h; delete hz;
   }

   auto corr = [](const std::vector<double> &X, const std::vector<double> &Y) {
      const int n = X.size();
      double sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0;
      for (int i = 0; i < n; ++i) { sx += X[i]; sy += Y[i]; sxx += X[i] * X[i];
         sxy += X[i] * Y[i]; syy += Y[i] * Y[i]; }
      return (n * sxy - sx * sy) / std::sqrt((n * sxx - sx * sx) * (n * syy - sy * sy));
   };
   auto spread = [](const std::vector<double> &V) {
      double mn = 1e9, mx = -1e9, s = 0;
      for (double v : V) { mn = std::min(mn, v); mx = std::max(mx, v); s += v; }
      const double mean = s / V.size();
      double sd = 0;
      for (double v : V) sd += (v - mean) * (v - mean);
      return std::make_tuple(mean, std::sqrt(sd / V.size()), mn, mx);
   };

   auto [mP, sP, nP, xP] = spread(exPk);
   auto [mM, sM, nM, xM] = spread(exMed);
   printf("\n=== matFX-OFF energy scale across %zu runs (unbiased by the collapse) ===\n", exPk.size());
   printf("  Ex peak   : mean %+.3f  rms %.3f  range [%+.3f, %+.3f] MeV\n", mP, sP, nP, xP);
   printf("  Ex median : mean %+.3f  rms %.3f  range [%+.3f, %+.3f] MeV\n", mM, sM, nM, xM);
   printf("\n=== does that scale explain the CATIMA collapse? ===\n");
   printf("  r(collapse, Ex peak)   = %+.3f\n", corr(col, exPk));
   printf("  r(collapse, Ex median) = %+.3f\n", corr(col, exMed));
   printf("  r(collapse, vertex z)  = %+.3f   [for reference]\n", corr(col, vz));
   printf("\nA flat Ex scale with a near-zero correlation kills the drift-velocity hypothesis:\n"
          "the runs that collapse are not geometrically different from the ones that do not.\n\n");

   auto *c = new TCanvas("c", "", 1500, 560);
   c->Divide(2, 1);
   c->cd(1);
   auto *g1 = new TGraph(col.size(), exPk.data(), col.data());
   g1->SetTitle("per-run geometry vs collapse;E_{x} peak, matFX-off sample  (MeV);CATIMA fits collapsed  (%)");
   g1->SetMarkerStyle(20); g1->SetMarkerSize(1.3); g1->SetMarkerColor(kAzure + 2);
   g1->Draw("AP");
   c->cd(2);
   auto *g2 = new TGraph(col.size(), vz.data(), col.data());
   g2->SetTitle("vertex z vs collapse;mean converged vertex z  (mm);CATIMA fits collapsed  (%)");
   g2->SetMarkerStyle(20); g2->SetMarkerSize(1.3); g2->SetMarkerColor(kRed + 1);
   g2->Draw("AP");
   c->SaveAs("plots/per_run_geometry.png");
   printf("wrote plots/per_run_geometry.png\n\n");
}
