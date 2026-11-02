/// @file order_score5_dp.C
/// @brief Repairs the data referee of order_score*_dp.C (2026-11-02): a fit collapsing to KE ~ 0 has Ex ~ 0 in elastic kinematics at any angle and
///        passes |Ex| < 2 MeV. GOOD now also requires the fitted Brho to agree with the fit-independent circle Brho, |ln(Brho_fit/bA)| < 0.5
///        (the circle is low by at most ~25 % for deuterons on sim). Reports how many previously GOOD tracks fail it.
///   root -l -b -q 'order_score5_dp.C("data0,data1",36860)'
namespace ords5 {
double exEl(double ke, double th, double tb)
{
   const double mC = 16.0147 * 931.49410, md = 1875.613;
   const double Eb = tb + mC, pb = std::sqrt(Eb * Eb - mC * mC), Ed = ke + md, pd = std::sqrt(Ed * Ed - md * md);
   const double E4 = Eb + md - Ed, px = -pd * std::sin(th), pz = pb - pd * std::cos(th);
   return std::sqrt(E4 * E4 - px * px - pz * pz) - mC;
}
double q(std::vector<double> v, double f)
{
   if (v.empty()) return -99;
   std::sort(v.begin(), v.end());
   return v[std::min(v.size() - 1, (size_t)(f * v.size()))];
}
}
void order_score5_dp(TString tags, Long64_t rowOff, TString O = "/mnt/f/a1975/order_test/")
{
   TChain a("o"), e("o"), g("o");
   for (auto *x : *tags.Tokenize(",")) { TString t = ((TObjString *)x)->GetString(); a.Add(O + t + ".root"); e.Add(O + "e_" + t + ".root"); g.Add(O + "g_" + t + ".root"); }
   const int NM = 3;
   const char *nm[NM] = {"0 production", "8 stray hits dropped", "9 dropped, large-pad guard"};
   int ok[NM]; double ke[NM], th[NM];
   a.SetBranchAddress("ok0", &ok[0]); a.SetBranchAddress("ke0", &ke[0]); a.SetBranchAddress("th0", &th[0]);
   e.SetBranchAddress("ok8", &ok[1]); e.SetBranchAddress("ke8", &ke[1]); e.SetBranchAddress("th8", &th[1]);
   g.SetBranchAddress("ok9", &ok[2]); g.SetBranchAddress("ke9", &ke[2]); g.SetBranchAddress("th9", &th[2]);
   TFile fp("/mnt/f/a1975/pid3d_sbpid/pv_gainnorm.root");
   auto *w = (TTree *)fp.Get("w");
   TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30&&s1_g>0&&bA>0";
   w->Draw("s1_g", sel, "goff");
   std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
   std::sort(x.begin(), x.end());
   const double md = x[x.size() / 2];
   const int NS = 8;
   const double bl[NS + 1] = {0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.45, 0.60, 0.80};
   double dm[NS];
   for (int k = 0; k < NS; ++k) {
      TH1D h("h", "", 600, -4, 2);
      w->Draw(Form("log(s1_g/%g*bA)>>h", md), sel + Form("&&bA>%g&&bA<%g", bl[k], bl[k + 1]), "goff");
      h.Smooth(2);
      dm[k] = h.GetBinCenter(h.GetMaximumBin());
   }
   double bA, s1g, sp3; int v3, ok1;
   w->SetBranchAddress("bA", &bA); w->SetBranchAddress("s1_g", &s1g); w->SetBranchAddress("sp3", &sp3); w->SetBranchAddress("v3", &v3); w->SetBranchAddress("ok1", &ok1);
   long nAll[NS] = {0}, gOld[NS][NM] = {{0}}, gNew[NS][NM] = {{0}}, tot = 0, tOld[NM] = {0}, tNew[NM] = {0};
   std::vector<double> E[NS][NM];
   for (Long64_t i = 0; i < a.GetEntries(); ++i) {
      a.GetEntry(i); e.GetEntry(i); g.GetEntry(i); w->GetEntry(i + rowOff);
      if (!v3 || !ok1 || sp3 >= 6 || s1g <= 0 || bA <= 0) continue;
      int k = -1;
      for (int j = 0; j < NS; ++j) if (bA > bl[j] && bA < bl[j + 1]) k = j;
      if (k < 0 || std::abs(std::log(s1g / md * bA) - dm[k]) > 0.12) continue;
      ++nAll[k]; ++tot;
      for (int m = 0; m < NM; ++m) {
         if (!ok[m]) continue;
         const double ex = ords5::exEl(ke[m], th[m], 180.0), b = std::sqrt(ke[m] * ke[m] + 2 * ke[m] * 1875.613) / 299.792458;
         const bool agree = std::abs(std::log(b / bA)) < 0.5;
         if (std::abs(ex) < 2) { ++gOld[k][m]; ++tOld[m]; if (agree) { ++gNew[k][m]; ++tNew[m]; } }
         if (agree) E[k][m].push_back(ex);
      }
   }
   for (int k = 0; k < NS; ++k) {
      if (nAll[k] < 15) continue;
      printf(" circle Brho %.2f-%.2f  N %ld\n", bl[k], bl[k + 1], nAll[k]);
      for (int m = 0; m < NM; ++m)
         printf("   %-28s GOOD old %5.1f %%  -> with Brho agreement %5.1f %%   Ex %+.2f [%.2f]\n", nm[m], 100. * gOld[k][m] / nAll[k], 100. * gNew[k][m] / nAll[k], ords5::q(E[k][m], .5), .5 * (ords5::q(E[k][m], .84) - ords5::q(E[k][m], .16)));
   }
   printf(" ALL %ld\n", tot);
   for (int m = 0; m < NM; ++m) printf("   %-28s GOOD old %5.1f %%  -> with Brho agreement %5.1f %%\n", nm[m], 100. * tOld[m] / tot, 100. * tNew[m] / tot);
}
