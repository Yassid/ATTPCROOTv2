/// @file order_score2_dp.C
/// @brief Separates what the re-clustering of order_fit_dp.C buys (2026-11-02). Modes 0 and 4 from the first pass, 5-7 (controls)
///        from the second, same tracks row by row. Sim: GOOD = fitted with |ln(Brho/true)| < 0.2, over ALL tracks. Data: tracks of the
///        deuteron band of the circle PID, GOOD = fitted with |Ex| < 2 MeV in elastic kinematics, over ALL of them.
///   root -l -b -q 'order_score2_dp.C("d5001,d5002,d5003,d5004",2.01355321,-1)'     sim
///   root -l -b -q 'order_score2_dp.C("data0,data1",2.01355321,36860)'               data (rowOff >= 0)
namespace ords2 {
double q(std::vector<double> v, double f)
{
   if (v.empty()) return -99;
   std::sort(v.begin(), v.end());
   return v[std::min(v.size() - 1, (size_t)(f * v.size()))];
}
double exEl(double ke, double th, double tb)
{
   const double mC = 16.0147 * 931.49410, md = 1875.613;
   const double Eb = tb + mC, pb = std::sqrt(Eb * Eb - mC * mC), Ed = ke + md, pd = std::sqrt(Ed * Ed - md * md);
   const double E4 = Eb + md - Ed, px = -pd * std::sin(th), pz = pb - pd * std::cos(th);
   return std::sqrt(E4 * E4 - px * px - pz * pz) - mC;
}
}
void order_score2_dp(TString tags, double massAmu, Long64_t rowOff, TString O = "/mnt/f/a1975/order_test/")
{
   const double mass = massAmu * 931.49410;
   TChain a("o"), b("o");
   for (auto *x : *tags.Tokenize(",")) { TString t = ((TObjString *)x)->GetString(); a.Add(O + t + ".root"); b.Add(O + "c_" + t + ".root"); }
   const int NM = 5;
   const int mode[NM] = {0, 4, 5, 6, 7};
   const char *nm[NM] = {"0 production", "4 re-clustered, 3D direction", "5 re-clustered, GeoTheta", "6 re-clustered, then z sort", "7 as 5, keeping blob hits"};
   int ok[NM], nh, nh2, skip;
   double ke[NM], th[NM], vr[NM], tbr;
   a.SetBranchAddress("nh", &nh); b.SetBranchAddress("nh", &nh2); a.SetBranchAddress("tbr", &tbr); a.SetBranchAddress("skip", &skip);
   for (int i = 0; i < NM; ++i) {
      TChain *c = mode[i] <= 4 ? &a : &b;
      c->SetBranchAddress(Form("ok%d", mode[i]), &ok[i]); c->SetBranchAddress(Form("ke%d", mode[i]), &ke[i]);
      c->SetBranchAddress(Form("th%d", mode[i]), &th[i]); c->SetBranchAddress(Form("vr%d", mode[i]), &vr[i]);
   }
   const bool isData = rowOff >= 0;
   TFile *fp = isData ? TFile::Open("/mnt/f/a1975/pid3d_sbpid/pv_gainnorm.root") : nullptr;
   TTree *w = isData ? (TTree *)fp->Get("w") : nullptr;
   const int NS = 8;
   const double bl[NS + 1] = {0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.45, 0.60, 0.80};
   double dm[NS] = {0}, md = 1, bA = 0, s1g = 0, sp3 = 0;
   int v3 = 0, ok1 = 0;
   if (isData) {
      TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30&&s1_g>0&&bA>0";
      w->Draw("s1_g", sel, "goff");
      std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
      std::sort(x.begin(), x.end());
      md = x[x.size() / 2];
      for (int k = 0; k < NS; ++k) {
         TH1D h("h", "", 600, -4, 2);
         w->Draw(Form("log(s1_g/%g*bA)>>h", md), sel + Form("&&bA>%g&&bA<%g", bl[k], bl[k + 1]), "goff");
         h.Smooth(2);
         dm[k] = h.GetBinCenter(h.GetMaximumBin());
      }
      w->SetBranchAddress("bA", &bA); w->SetBranchAddress("s1_g", &s1g); w->SetBranchAddress("sp3", &sp3); w->SetBranchAddress("v3", &v3); w->SetBranchAddress("ok1", &ok1);
   }
   std::vector<double> E[NS][NM];
   long nAll[NS] = {0}, nFit[NS][NM] = {{0}}, nGood[NS][NM] = {{0}}, mism = 0;
   for (Long64_t i = 0; i < a.GetEntries(); ++i) {
      a.GetEntry(i); b.GetEntry(i);
      if (nh != nh2) { ++mism; continue; }
      int k = -1;
      if (isData) {
         w->GetEntry(i + rowOff);
         if (!v3 || !ok1 || sp3 >= 6 || s1g <= 0 || bA <= 0) continue;
         for (int j = 0; j < NS; ++j) if (bA > bl[j] && bA < bl[j + 1]) k = j;
         if (k < 0 || std::abs(std::log(s1g / md * bA) - dm[k]) > 0.12) continue;
      } else {
         if (skip) continue;
         for (int j = 0; j < NS; ++j) if (tbr > bl[j] && tbr < bl[j + 1]) k = j;
         if (k < 0) continue;
      }
      ++nAll[k];
      for (int m = 0; m < NM; ++m) {
         if (!ok[m]) continue;
         ++nFit[k][m];
         double v = isData ? ords2::exEl(ke[m], th[m], 180.0) : std::log(std::sqrt(ke[m] * ke[m] + 2 * ke[m] * mass) / 299.792458 / tbr);
         E[k][m].push_back(v);
         nGood[k][m] += std::abs(v) < (isData ? 2.0 : 0.2);
      }
   }
   printf("mismatches %ld; %s\n", mism, isData ? "Ex [MeV], GOOD |Ex|<2" : "ln(Brho/true), GOOD |.|<0.2");
   long tot = 0, tG[NM] = {0}, tF[NM] = {0};
   for (int k = 0; k < NS; ++k) {
      if (nAll[k] < 15) continue;
      printf(" Brho %.2f-%.2f  N %ld\n", bl[k], bl[k + 1], nAll[k]);
      tot += nAll[k];
      for (int m = 0; m < NM; ++m) {
         tG[m] += nGood[k][m]; tF[m] += nFit[k][m];
         printf("   %-30s fit %5.1f %%  median %+.3f [%.3f]  GOOD %5.1f %%\n", nm[m], 100. * nFit[k][m] / nAll[k], ords2::q(E[k][m], .5), .5 * (ords2::q(E[k][m], .84) - ords2::q(E[k][m], .16)), 100. * nGood[k][m] / nAll[k]);
      }
   }
   printf(" ALL N %ld\n", tot);
   for (int m = 0; m < NM; ++m) printf("   %-30s fit %5.1f %%  GOOD %5.1f %%\n", nm[m], 100. * tF[m] / tot, 100. * tG[m] / tot);
}
