/// @file ab3d_physics_dp.C
/// @brief Physics-level score of the (d,p) direction A/B (default AtPRA direction "c" vs AtHelix3D "s"),
///        with the ADOPTED explorer selection applied through exsel_dp.C's own functions (chi2, angle-
///        corrected theta window, KE+keOff, vz window, two-body Ex). Input: dir_ab_dp.C's tree "ab".
///        Reports, per arm: selected tracks, bound-region counts (Ex -1.5..1.0) by theta_cm, the g.s.
///        peak (gaussian + linear, Poisson likelihood) and the same restricted to the tracks whose fit
///        the 3D direction CHANGED. Writes an Ex overlay PNG.
///   root -l -b -q 'ab3d_physics_dp.C("/mnt/f/a1975/gf_dp_ab3d/ab.root","plots/ab3d_physics_dp.png")'
#include "exsel_dp.C"

void ab3d_physics_dp(TString in, TString png, double vzLo = -100, double vzHi = 940)
{
   using namespace exsel;
   gStyle->SetOptStat(0);
   Cuts S;
   S.vzLo = vzLo; S.vzHi = vzHi;
   TFile f(in);
   auto *t = (TTree *)f.Get("ab");
   int F[2];
   double KE[2], TH[2], C2[2], VZ[2];
   const char *arm[2] = {"c", "s"}, *nm[2] = {"AtPRA direction (default)", "AtHelix3D direction"};
   for (int a = 0; a < 2; ++a) {
      t->SetBranchAddress(Form("%sfit", arm[a]), &F[a]); t->SetBranchAddress(Form("%ske", arm[a]), &KE[a]);
      t->SetBranchAddress(Form("%sth", arm[a]), &TH[a]); t->SetBranchAddress(Form("%schi", arm[a]), &C2[a]);
      t->SetBranchAddress(Form("%svz", arm[a]), &VZ[a]);
   }
   TH1D *hEx[2], *hExCh[2];
   const int tb[5] = {10, 20, 30, 40, 50};
   long nSel[2] = {0, 0}, nBound[2][4] = {{0}}, nBoundAll[2] = {0, 0}, nCh = 0;
   for (int a = 0; a < 2; ++a) {
      hEx[a] = new TH1D(Form("hEx%d", a), "16C(d,p)17C, adopted selection;E_{x} [MeV];counts / 155 keV", S.exBins, S.exLo, S.exHi);
      hExCh[a] = new TH1D(Form("hExCh%d", a), "only tracks whose fit CHANGED;E_{x} [MeV];counts", 45, S.exLo, S.exHi);
   }
   auto select = [&](int a, double &ex, double &tcm) {
      if (!F[a]) return false;
      const double k = KE[a] + S.keOff, thc = thCorr(S, TH[a], k);
      if (!(C2[a] <= S.chi2) || !(thc >= S.thLo && thc <= S.thHi) || !(k >= S.keLo && k <= S.keHi) ||
          !(VZ[a] >= S.vzLo && VZ[a] <= S.vzHi))
         return false;
      return kine2b(S.ebeam, thc * TMath::Pi() / 180.0, k, ex, tcm);
   };
   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      bool changed = F[0] != F[1] || (F[0] && F[1] && (std::abs(KE[0] - KE[1]) > 1e-6 || std::abs(TH[0] - TH[1]) > 1e-6));
      nCh += changed;
      for (int a = 0; a < 2; ++a) {
         double ex, tcm;
         if (!select(a, ex, tcm)) continue;
         ++nSel[a];
         hEx[a]->Fill(ex);
         if (changed) hExCh[a]->Fill(ex);
         if (ex > -1.5 && ex < 1.0) {
            ++nBoundAll[a];
            for (int b = 0; b < 4; ++b)
               if (tcm >= tb[b] && tcm < tb[b + 1]) ++nBound[a][b];
         }
      }
   }
   printf("tracks in the A/B tree %lld; fit changed by the 3D direction: %ld (%.1f %%)\n", t->GetEntries(), nCh,
          100.0 * nCh / t->GetEntries());
   printf("adopted selection, vz %g..%g:\n", vzLo, vzHi);
   printf("  %-28s %8s %8s | bound (Ex -1.5..1.0) by theta_cm: 10-20 20-30 30-40 40-50\n", "", "selected", "bound");
   for (int a = 0; a < 2; ++a)
      printf("  %-28s %8ld %8ld |                              %5ld %5ld %5ld %5ld\n", nm[a], nSel[a], nBoundAll[a],
             nBound[a][0], nBound[a][1], nBound[a][2], nBound[a][3]);
   // g.s. peak: gaussian + linear background, Poisson likelihood, Ex -3..2
   for (int a = 0; a < 2; ++a) {
      TF1 g(Form("g%d", a), "gaus(0)+pol1(3)", -3, 2);
      g.SetParameters(hEx[a]->GetMaximum(), -0.3, 0.3, 1, 0);
      g.SetParLimits(2, 0.08, 1.0);
      auto r = hEx[a]->Fit(&g, "QNRLS");
      double area = g.GetParameter(0) * g.GetParameter(2) * std::sqrt(2 * M_PI) / hEx[a]->GetBinWidth(1);
      printf("  g.s. fit %-28s mean %+.3f  sigma %.3f +- %.3f MeV  peak counts %.0f  (status %d)\n", nm[a],
             g.GetParameter(1), g.GetParameter(2), g.GetParError(2), area, (int)r);
   }
   printf("  changed tracks only: selected %0.f (default) vs %0.f (3D); bound %0.f vs %0.f\n", hExCh[0]->Integral(),
          hExCh[1]->Integral(), hExCh[0]->Integral(hExCh[0]->FindBin(-1.49), hExCh[0]->FindBin(0.99)),
          hExCh[1]->Integral(hExCh[1]->FindBin(-1.49), hExCh[1]->FindBin(0.99)));
   auto *c = new TCanvas("c", "c", 1400, 600);
   c->Divide(2, 1);
   int col[2] = {kBlack, kRed + 1};
   for (int p = 0; p < 2; ++p) {
      c->cd(p + 1);
      TH1D **h = p == 0 ? hEx : hExCh;
      double mx = std::max(h[0]->GetMaximum(), h[1]->GetMaximum()) * 1.15;
      for (int a = 0; a < 2; ++a) {
         h[a]->SetLineColor(col[a]); h[a]->SetLineWidth(2); h[a]->SetMaximum(mx);
         h[a]->Draw(a ? "hist same" : "hist");
      }
      auto *lg = new TLegend(0.55, 0.75, 0.89, 0.89);
      for (int a = 0; a < 2; ++a) lg->AddEntry(h[a], nm[a], "l");
      lg->Draw();
   }
   c->SaveAs(png);
}
