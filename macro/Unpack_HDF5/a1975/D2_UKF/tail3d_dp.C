/// @file tail3d_dp.C
/// @brief Referee deuterons (pid_dd_referee.C sample) that the 3D helix-path PID puts in the PROTON band
///        while Spyral does not: which ingredient is off? Splits ln m_3D - ln m_Spyral into Brho and
///        dE/dx parts and compares the 3D measurement's internals (radius, helix path, Spyral arc,
///        window) between this tail and the well-measured deuterons.
///   root -l -b -q 'tail3d_dp.C("/mnt/f/a1975/reco_d2_tau/pid_width_0016.root")'
#include "pid_dd_referee.C"
void tail3d_dp(TString in, double Eb = 184.25)
{
   TFile f(in);
   auto *w = (TTree *)f.Get("w");
   TString S = "va&&v3&&sp3<6&&nh>=30&&sq3h>0";
   auto med = [&](const char *v) {
      w->Draw(v, S, "goff");
      std::vector<double> x(w->GetV1(), w->GetV1() + w->GetSelectedRows());
      std::sort(x.begin(), x.end());
      return x[x.size() / 2];
   };
   double mS = med("sqa"), mH = med("sq3h");
   int va, v3, nh, how;
   Long64_t ev;
   int tid;
   double br, b3, sqa, sq3h, pol, th3, sp3, pathH, pathF, spyArc, sWin3, R3, dca3;
   w->SetBranchAddress("ev", &ev); w->SetBranchAddress("tid", &tid);
   w->SetBranchAddress("va", &va); w->SetBranchAddress("v3", &v3); w->SetBranchAddress("nh", &nh);
   w->SetBranchAddress("how3", &how); w->SetBranchAddress("br", &br); w->SetBranchAddress("b3", &b3);
   w->SetBranchAddress("sqa", &sqa); w->SetBranchAddress("sq3h", &sq3h); w->SetBranchAddress("pol", &pol);
   w->SetBranchAddress("th3", &th3); w->SetBranchAddress("sp3", &sp3); w->SetBranchAddress("pathH", &pathH);
   w->SetBranchAddress("pathF", &pathF); w->SetBranchAddress("spyArc", &spyArc); w->SetBranchAddress("sWin3", &sWin3);
   w->SetBranchAddress("R3", &R3); w->SetBranchAddress("dca3", &dca3);
   enum { kB, kE, kPathRatio, kArcRatio, kTh, kNh, kR, kDca, kPathH, kN };
   const char *nm[kN] = {"ln(Brho 3D/Spyral)", "ln(sqrt dEdx 3D/Spyral)", "helix path / centroid path",
                         "helix path / Spyral arc", "theta_lab 3D", "hits", "R vertex-end helix", "dca",
                         "helix path [mm]"};
   std::vector<double> T[kN], G[kN];
   std::vector<std::pair<Long64_t, int>> list;
   for (Long64_t i = 0; i < w->GetEntries(); ++i) {
      w->GetEntry(i);
      if (!va || !v3 || sp3 >= 6 || nh < 30 || sq3h <= 0) continue;
      if (!(std::abs(exDD(Eb, (180 - pol) * M_PI / 180, br)) < 1.5 && std::abs(exDD(Eb, th3 * M_PI / 180, b3)) < 1.5))
         continue;
      double lS = std::log(sqa / mS * br), lH = std::log(sq3h / mH * b3);
      bool tH = lH < -1.470 - 0.45, tS = lS < -1.429 - 0.45;
      if (tS) continue;          // keep only deuterons Spyral measures as deuterons
      auto &V = tH ? T : G;
      V[kB].push_back(std::log(b3 / br));
      V[kE].push_back(std::log((sq3h / mH) / (sqa / mS)));
      V[kPathRatio].push_back(pathF > 0 ? pathH / pathF : -1);
      V[kArcRatio].push_back(spyArc > 0 ? pathH / spyArc : -1);
      V[kTh].push_back(th3); V[kNh].push_back(nh); V[kR].push_back(R3); V[kDca].push_back(dca3);
      V[kPathH].push_back(pathH);
      if (tH) list.push_back({ev, tid});
   }
   printf("deuterons Spyral calls deuterons: %zu; of them 3D puts %zu in the proton band\n", T[0].size() + G[0].size(),
          T[0].size());
   printf("  %-30s %-30s %-30s\n", "", "3D-ONLY TAIL  p16 / med / p84", "WELL MEASURED  p16 / med / p84");
   for (int k = 0; k < kN; ++k) {
      auto q = [](std::vector<double> x) {
         std::sort(x.begin(), x.end());
         int n = x.size();
         return n > 3 ? std::array<double, 3>{x[n * 16 / 100], x[n / 2], x[n * 84 / 100]} : std::array<double, 3>{0, 0, 0};
      };
      auto a = q(T[k]), b = q(G[k]);
      printf("  %-30s %8.3f %8.3f %8.3f       %8.3f %8.3f %8.3f\n", nm[k], a[0], a[1], a[2], b[0], b[1], b[2]);
   }
   printf("first tail tracks (ev:tid): ");
   for (size_t i = 0; i < list.size(); ++i) printf("%lld:%d ", list[i].first, list[i].second);
   printf("\n");
}
