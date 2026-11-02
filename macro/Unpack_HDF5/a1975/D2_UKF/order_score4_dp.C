/// @file order_score4_dp.C
/// @brief Large-pad guard (2026-11-02): modes 0 (first pass), 8 (third pass) and 9 (fourth pass) of order_fit_dp.C, same tracks row
///        by row, split by whether AtHelix3D dropped hits (nDrop). Sim: GOOD = fitted with |ln(Brho/true)| < 0.2, over ALL tracks. Data: tracks of the
///        deuteron band of the circle PID, GOOD = fitted with |Ex| < 2 MeV in elastic kinematics, over ALL of them.
///   root -l -b -q 'order_score4_dp.C("d5001,d5002,d5003,d5004",2.01355321,-1)'     sim
///   root -l -b -q 'order_score4_dp.C("data0,data1",2.01355321,36860)'               data (rowOff >= 0)
namespace ords4 {
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
void order_score4_dp(TString tags, double massAmu, Long64_t rowOff, TString O = "/mnt/f/a1975/order_test/")
{
   const double mass = massAmu * 931.49410;
   TChain a("o"), b("o"), e("o"), g("o");
   for (auto *x : *tags.Tokenize(",")) { TString t = ((TObjString *)x)->GetString(); a.Add(O + t + ".root"); b.Add(O + "c_" + t + ".root"); e.Add(O + "e_" + t + ".root"); g.Add(O + "g_" + t + ".root"); }
   const int NM = 3;
   const int mode[NM] = {0, 8, 9};
   const char *nm[NM] = {"0 production", "8 stray hits dropped", "9 dropped, large-pad guard"};
   int nDrop;
   e.SetBranchAddress("nDrop", &nDrop);
   int ok[NM], nh, nh2, skip;
   double ke[NM], th[NM], vr[NM], tbr;
   a.SetBranchAddress("nh", &nh); b.SetBranchAddress("nh", &nh2); a.SetBranchAddress("tbr", &tbr); a.SetBranchAddress("skip", &skip);
   for (int i = 0; i < NM; ++i) {
      TChain *c = mode[i] == 0 ? &a : (mode[i] == 8 ? &e : &g);
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
   std::vector<double> E[NS][NM], fd[NS];
   long nAll[NS] = {0}, nFit[NS][NM] = {{0}}, nGood[NS][NM] = {{0}}, mism = 0, nD[NS] = {0}, gD[NS][NM] = {{0}}, gK[NS][NM] = {{0}};
   for (Long64_t i = 0; i < a.GetEntries(); ++i) {
      a.GetEntry(i); b.GetEntry(i); e.GetEntry(i); g.GetEntry(i);
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
      const bool dropped = nDrop > 0;
      nD[k] += dropped;
      if (dropped) fd[k].push_back(double(nDrop) / nh);
      for (int m = 0; m < NM; ++m) {
         if (!ok[m]) continue;
         ++nFit[k][m];
         double v = isData ? ords4::exEl(ke[m], th[m], 180.0) : std::log(std::sqrt(ke[m] * ke[m] + 2 * ke[m] * mass) / 299.792458 / tbr);
         E[k][m].push_back(v);
         const bool g = std::abs(v) < (isData ? 2.0 : 0.2);
         nGood[k][m] += g;
         (dropped ? gD : gK)[k][m] += g;
      }
   }
   printf("mismatches %ld; %s\n", mism, isData ? "Ex [MeV], GOOD |Ex|<2" : "ln(Brho/true), GOOD |.|<0.2");
   long tot = 0, tG[NM] = {0}, tF[NM] = {0};
   for (int k = 0; k < NS; ++k) {
      if (nAll[k] < 15) continue;
      printf(" Brho %.2f-%.2f  N %ld; tracks with dropped hits %.0f %% (median share dropped %.0f %%)\n", bl[k], bl[k + 1], nAll[k], 100. * nD[k] / nAll[k], 100 * ords4::q(fd[k], .5));
      tot += nAll[k];
      for (int m = 0; m < NM; ++m) {
         tG[m] += nGood[k][m]; tF[m] += nFit[k][m];
         printf("   %-30s fit %5.1f %%  median %+.3f [%.3f]  GOOD %5.1f %%  (dropped: %5.1f %%, none dropped: %5.1f %%)\n", nm[m], 100. * nFit[k][m] / nAll[k], ords4::q(E[k][m], .5), .5 * (ords4::q(E[k][m], .84) - ords4::q(E[k][m], .16)), 100. * nGood[k][m] / nAll[k], 100. * gD[k][m] / std::max(1L, nD[k]), 100. * gK[k][m] / std::max(1L, nAll[k] - nD[k]));
      }
   }
   printf(" ALL N %ld\n", tot);
   for (int m = 0; m < NM; ++m) printf("   %-30s fit %5.1f %%  GOOD %5.1f %%\n", nm[m], 100. * tF[m] / tot, 100. * tG[m] / tot);
}
