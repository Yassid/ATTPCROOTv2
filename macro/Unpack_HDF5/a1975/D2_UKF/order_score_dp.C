/// @file order_score_dp.C
/// @brief Scores order_fit_dp.C against truth: per true-Brho slice and for near-90-degree tracks, for each cluster
///        ordering: fit rate, median ln(Brho/true) [half 68 %], gross failures |ln(Brho/true)| > 0.2, median chi2/ndf, polar
///        error. Only tracks fitted by the production mode 0 enter the comparison of fits ("common"), so a mode cannot look
///        better by losing its hard tracks; the fit rate is reported separately. Also the z-sort path inflation.
///   root -l -b -q 'order_score_dp.C("/mnt/f/a1975/order_test/d5001.root,...",2.01355321)'
namespace ordscore {
double q(std::vector<double> v, double f)
{
   if (v.empty()) return -99;
   std::sort(v.begin(), v.end());
   return v[std::min(v.size() - 1, (size_t)(f * v.size()))];
}
}
void order_score_dp(TString files, double massAmu)
{
   const double mass = massAmu * 931.49410;
   TChain o("o");
   for (auto *x : *files.Tokenize(",")) o.Add(((TObjString *)x)->GetString());
   const int NM = 5;
   const char *nm[NM] = {"0 z sort (production)", "1 z sort, 3D direction", "2 stored order", "3 3D backbone", "4 backbone, re-clustered"};
   double tbr, tth, pz, ps, pb, ke[NM], th[NM], c2[NM], nd[NM], vr[NM];
   int nh, ncl, skip, ok[NM];
   o.SetBranchAddress("tbr", &tbr); o.SetBranchAddress("tth", &tth); o.SetBranchAddress("pz", &pz); o.SetBranchAddress("ps", &ps);
   o.SetBranchAddress("pb", &pb); o.SetBranchAddress("nh", &nh); o.SetBranchAddress("ncl", &ncl); o.SetBranchAddress("skip", &skip);
   for (int m = 0; m < NM; ++m) {
      o.SetBranchAddress(Form("ok%d", m), &ok[m]); o.SetBranchAddress(Form("ke%d", m), &ke[m]); o.SetBranchAddress(Form("th%d", m), &th[m]);
      o.SetBranchAddress(Form("c2%d", m), &c2[m]); o.SetBranchAddress(Form("nd%d", m), &nd[m]); o.SetBranchAddress(Form("vr%d", m), &vr[m]);
   }
   // reco polar vs truth: decide the handedness once from mode 0
   double same = 0, flip = 0;
   for (Long64_t i = 0; i < o.GetEntries(); ++i) { o.GetEntry(i); if (skip || !ok[0]) continue; same += std::abs(th[0] - tth) < std::abs(th[0] - (M_PI - tth)); flip += std::abs(th[0] - tth) >= std::abs(th[0] - (M_PI - tth)); }
   const bool flipT = flip > same;
   printf("%lld rows; reco polar compares with %s\n", o.GetEntries(), flipT ? "180-true" : "true");
   const int NS = 8;
   const double bl[NS + 1] = {0.10, 0.20, 0.25, 0.30, 0.35, 0.45, 0.60, 0.80, 1.4};
   for (int sel = 0; sel < 2; ++sel) {
      printf("\n=== %s\n", sel == 0 ? "ALL POLAR ANGLES" : "NEAR 90 deg: |theta_true - 90| < 5");
      for (int k = 0; k < NS; ++k) {
         std::vector<double> B[NM], C[NM], T[NM], infl;
         long nAll = 0, nFit[NM] = {0};
         for (Long64_t i = 0; i < o.GetEntries(); ++i) {
            o.GetEntry(i);
            if (skip || tbr < bl[k] || tbr >= bl[k + 1]) continue;
            const double tdeg = tth * 57.29578;
            if (sel == 1 && std::abs(tdeg - 90) > 5) continue;
            ++nAll;
            if (ps > 0 && pz > 0) infl.push_back(pz / ps);
            for (int m = 0; m < NM; ++m) nFit[m] += ok[m];
            if (!ok[0]) continue;
            for (int m = 0; m < NM; ++m) {
               if (!ok[m]) continue;
               const double br = std::sqrt(ke[m] * ke[m] + 2 * ke[m] * mass) / 299.792458;
               B[m].push_back(std::log(br / tbr));
               C[m].push_back(c2[m] / nd[m]);
               T[m].push_back((th[m] - (flipT ? M_PI - tth : tth)) * 57.29578);
            }
         }
         if (nAll < 15) continue;
         printf(" true Brho %.2f-%.2f  N %ld  z-sort path / stored path: median %.2f p84 %.2f, >1.5: %.0f %%\n", bl[k], bl[k + 1], nAll, ordscore::q(infl, .5), ordscore::q(infl, .84),
                100. * std::count_if(infl.begin(), infl.end(), [](double x) { return x > 1.5; }) / std::max<size_t>(1, infl.size()));
         for (int m = 0; m < NM; ++m) {
            long gross = 0;
            for (double b : B[m]) gross += std::abs(b) > 0.2;
            printf("   %-26s fit %5.1f %%  ln(B/true) %+.3f [%.3f]  gross %5.1f %%  chi2/ndf %7.3f  dtheta %+.2f [%.2f]\n", nm[m], 100. * nFit[m] / nAll, ordscore::q(B[m], .5),
                   .5 * (ordscore::q(B[m], .84) - ordscore::q(B[m], .16)), 100. * gross / std::max<size_t>(1, B[m].size()), ordscore::q(C[m], .5), ordscore::q(T[m], .5), .5 * (ordscore::q(T[m], .84) - ordscore::q(T[m], .16)));
         }
      }
   }
}
