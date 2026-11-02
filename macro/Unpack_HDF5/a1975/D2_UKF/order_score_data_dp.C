/// @file order_score_data_dp.C
/// @brief Scores order_fit_dp.C on DATA with a referee independent of the fit (2026-11-02): tracks in the deuteron band of the
///        circle PID (pv tree, |ln-mass| < 0.12, chosen before any fit), fitted as deuterons; elastic d(16C,d) kinematics require
///        Ex = 0 for every one of them. Per circle-Brho slice and mode: fit rate; median Ex [half 68 %] of the fitted; and
///        GOOD = fitted with |Ex| < exCut as a fraction of ALL band tracks, which a mode cannot improve by losing hard tracks.
///        Same statistics for vertex distance from the axis. Plus the z-sort path inflation.
///   root -l -b -q 'order_score_data_dp.C("/mnt/f/a1975/order_test/data0.root,/mnt/f/a1975/order_test/data1.root",36860)'
namespace ordsd {
double q(std::vector<double> v, double f)
{
   if (v.empty()) return -99;
   std::sort(v.begin(), v.end());
   return v[std::min(v.size() - 1, (size_t)(f * v.size()))];
}
/// Ex of 16C from a deuteron of kinetic energy ke at polar angle th, beam kinetic energy tb
double exEl(double ke, double th, double tb)
{
   const double mC = 16.0147 * 931.49410, md = 1875.613;
   const double Eb = tb + mC, pb = std::sqrt(Eb * Eb - mC * mC), Ed = ke + md, pd = std::sqrt(Ed * Ed - md * md);
   const double E4 = Eb + md - Ed, px = -pd * std::sin(th), pz = pb - pd * std::cos(th);
   return std::sqrt(E4 * E4 - px * px - pz * pz) - mC;
}
}
void order_score_data_dp(TString files, Long64_t rowOff, double exCut = 2.0, double tBeam = 180.0, TString pv = "/mnt/f/a1975/pid3d_sbpid/pv_gainnorm.root")
{
   TChain o("o");
   for (auto *x : *files.Tokenize(",")) o.Add(((TObjString *)x)->GetString());
   const int NM = 5;
   const char *nm[NM] = {"0 z sort (production)", "1 z sort, 3D direction", "2 stored order", "3 3D backbone", "4 backbone, re-clustered"};
   double pz, ps, ke[NM], th[NM], c2[NM], nd[NM], vr[NM];
   int nh, ok[NM];
   o.SetBranchAddress("pz", &pz); o.SetBranchAddress("ps", &ps); o.SetBranchAddress("nh", &nh);
   for (int m = 0; m < NM; ++m) {
      o.SetBranchAddress(Form("ok%d", m), &ok[m]); o.SetBranchAddress(Form("ke%d", m), &ke[m]); o.SetBranchAddress(Form("th%d", m), &th[m]);
      o.SetBranchAddress(Form("c2%d", m), &c2[m]); o.SetBranchAddress(Form("nd%d", m), &nd[m]); o.SetBranchAddress(Form("vr%d", m), &vr[m]);
   }
   TFile fp(pv);
   auto *w = (TTree *)fp.Get("w");
   TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30&&s1_g>0&&bA>0";
   w->Draw("s1_g", sel, "goff");
   std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
   std::sort(x.begin(), x.end());
   const double md = x[x.size() / 2];
   const int NS = 7;
   const double bl[NS + 1] = {0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.45, 0.60};
   double dm[NS];
   for (int k = 0; k < NS; ++k) {
      TH1D h("h", "", 600, -4, 2);
      w->Draw(Form("log(s1_g/%g*bA)>>h", md), sel + Form("&&bA>%g&&bA<%g", bl[k], bl[k + 1]), "goff");
      h.Smooth(2);
      dm[k] = h.GetBinCenter(h.GetMaximumBin());
   }
   double bA, s1g, sp3;
   int va, v3, ok1, nhp;
   w->SetBranchAddress("bA", &bA); w->SetBranchAddress("s1_g", &s1g); w->SetBranchAddress("sp3", &sp3); w->SetBranchAddress("va", &va);
   w->SetBranchAddress("v3", &v3); w->SetBranchAddress("ok1", &ok1); w->SetBranchAddress("nh", &nhp);
   std::vector<double> E[NS][NM], V[NS][NM], infl[NS];
   long nAll[NS] = {0}, nFit[NS][NM] = {{0}}, nGood[NS][NM] = {{0}}, mism = 0;
   for (Long64_t i = 0; i < o.GetEntries(); ++i) {
      o.GetEntry(i);
      w->GetEntry(i + rowOff);
      if (nhp != nh) { ++mism; continue; }
      if (!v3 || !ok1 || sp3 >= 6 || s1g <= 0 || bA <= 0) continue;
      int k = -1;
      for (int j = 0; j < NS; ++j) if (bA > bl[j] && bA < bl[j + 1]) k = j;
      if (k < 0 || std::abs(std::log(s1g / md * bA) - dm[k]) > 0.12) continue;
      ++nAll[k];
      if (ps > 0 && pz > 0) infl[k].push_back(pz / ps);
      for (int m = 0; m < NM; ++m) {
         if (!ok[m]) continue;
         ++nFit[k][m];
         const double ex = ordsd::exEl(ke[m], th[m], tBeam);
         E[k][m].push_back(ex);
         V[k][m].push_back(vr[m]);
         nGood[k][m] += std::abs(ex) < exCut;
      }
   }
   printf("rows %lld, row mismatches %ld; beam %.0f MeV, good = |Ex| < %.1f MeV\n", o.GetEntries(), mism, tBeam, exCut);
   long tot = 0, totG[NM] = {0}, totF[NM] = {0};
   for (int k = 0; k < NS; ++k) {
      if (nAll[k] < 15) continue;
      long big = 0;
      for (double v : infl[k]) big += v > 1.5;
      printf(" circle Brho %.2f-%.2f  deuteron-band tracks %ld; z-sort path/stored path median %.2f, >1.5: %.0f %%\n", bl[k], bl[k + 1], nAll[k], ordsd::q(infl[k], .5), 100. * big / std::max<size_t>(1, infl[k].size()));
      tot += nAll[k];
      for (int m = 0; m < NM; ++m) {
         totG[m] += nGood[k][m]; totF[m] += nFit[k][m];
         printf("   %-26s fit %5.1f %%  Ex %+6.2f [%5.2f]  GOOD %5.1f %%  vertex |xy| %6.1f\n", nm[m], 100. * nFit[k][m] / nAll[k], ordsd::q(E[k][m], .5), .5 * (ordsd::q(E[k][m], .84) - ordsd::q(E[k][m], .16)), 100. * nGood[k][m] / nAll[k], ordsd::q(V[k][m], .5));
      }
   }
   printf(" ALL %ld tracks:\n", tot);
   for (int m = 0; m < NM; ++m) printf("   %-26s fit %5.1f %%  GOOD %5.1f %%\n", nm[m], 100. * totF[m] / tot, 100. * totG[m] / tot);
}
