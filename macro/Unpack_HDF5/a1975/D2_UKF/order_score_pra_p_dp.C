/// @file order_score_pra_p_dp.C
/// @brief Backward (d,p) protons with AtPRA::SetRemoveStrayHits (2026-11-02): proton-band tracks of runs 0016 0026 0058 0088, same selection and Ex as
///        order_score_p_dp.C; production clusters vs the option, matched by (event, trackID). Compare with the find805 line of order_score_p_dp.C.
///   root -l -b -q 'order_score_pra_p_dp.C()'
namespace ordpp {
double exDP(double ke, double th, double tb)
{
   const double mC = 16.0147 * 931.49410, md = 1875.613, mp = 938.272, M4 = mC + md - mp + 1.496;
   const double Eb = tb + mC, pb = std::sqrt(Eb * Eb - mC * mC), Ep = ke + mp, pp = std::sqrt(Ep * Ep - mp * mp);
   const double E4 = Eb + md - Ep, px = -pp * std::sin(th), pz = pb - pp * std::cos(th);
   return std::sqrt(E4 * E4 - px * px - pz * pz) - M4;
}
}
void order_score_pra_p_dp(TString O = "/mnt/f/a1975/order_test/", TString Q = "/mnt/f/a1975/pra_test/")
{
   const char *runs[4] = {"0016", "0026", "0058", "0088"};
   long sel[2] = {0, 0}, bnd[2] = {0, 0}, junk[2] = {0, 0}, off[2] = {0, 0}, fit[2] = {0, 0}, all[2] = {0, 0};
   for (int r = 0; r < 4; ++r) {
      TChain a("o");
      a.Add(O + "p_" + runs[r] + "_0.root"); a.Add(O + "p_" + runs[r] + "_1.root");
      TChain b("o");
      b.Add(Q + "f_p" + runs[r] + ".root");
      TTree *T[2] = {&a, &b};
      for (int m = 0; m < 2; ++m) {
         int ok; double ke, th, c2, nd, vr, vz;
         T[m]->SetBranchAddress("ok0", &ok); T[m]->SetBranchAddress("ke0", &ke); T[m]->SetBranchAddress("th0", &th); T[m]->SetBranchAddress("c20", &c2);
         T[m]->SetBranchAddress("nd0", &nd); T[m]->SetBranchAddress("vr0", &vr); T[m]->SetBranchAddress("vz0", &vz);
         for (Long64_t i = 0; i < T[m]->GetEntries(); ++i) {
            T[m]->GetEntry(i);
            ++all[m];
            if (!ok) continue;
            ++fit[m];
            const double deg = th * 57.29578, chi = c2 / nd;
            if (deg < 89 || deg > 161 || ke <= 1.2 || chi >= 13.5 || vz < -100 || vz > 940) continue;
            const double ex = ordpp::exDP(ke + 0.25, th, 184.25 - 0.01103 * (vz - 500));
            ++sel[m];
            if (ex > -1.5 && ex < 1.0) ++bnd[m];
            if (ex < -2) ++junk[m];
            if (vr > 20) ++off[m];
         }
      }
   }
   const char *nm[2] = {"production clusters", "AtPRA::SetRemoveStrayHits"};
   for (int m = 0; m < 2; ++m)
      printf("  %-28s tracks %5ld fit %5.1f %%  selected %5ld  bound %4ld  junk %4ld (%.1f %%)  off-axis %4ld (%.1f %%)\n", nm[m], all[m], 100. * fit[m] / std::max(1L, all[m]), sel[m], bnd[m], junk[m], 100. * junk[m] / std::max(1L, sel[m]), off[m], 100. * off[m] / std::max(1L, sel[m]));
}
