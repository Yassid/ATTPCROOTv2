/// @file plot_pid_sim_methods_dp.C
/// @brief Simulated (d,t) tritons (g.s. + ex2, pid_sbpid_dp.C on sim) through every Brho / PID method (2026-11-01).
///        (a) Brho bias ln(B/B_true) and (b) its resolution vs TRUE Brho, for Spyral (br), 3D (b1), 3D hit-refit +
///        arc-length angle (bA), species-blind EOM combined (efc). (c) Where the tritons land relative to the DEUTERON
///        locus, delta = ln(hs / hs_d(B)) with hs_d the CATIMA deuteron prediction at that Brho (k_sim from truth),
///        vs each method's OWN reconstructed Brho (what a gate sees); the band is the data deuteron half-width (~0.08).
///        (d) the species-blind fitted mass ln(m_fit / m_t) vs true Brho (d would sit at -0.405).
///   root -l -b -q 'plot_pid_sim_methods_dp.C("plots/pid_sim_methods.png")'
#include "eom_brho_dp.h"
void plot_pid_sim_methods_dp(TString png, double ksim = 1.982e-5,
                             TString files = "/mnt/f/a1975/pid3d_sbpid/sim_gs.root,/mnt/f/a1975/pid3d_sbpid/sim_ex2.root")
{
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   eom::Table tab[3];
   eom::makeTables(tab);
   TChain w("w");
   for (auto *o : *files.Tokenize(",")) w.Add(((TObjString *)o)->GetString());
   double hs, bA, br, b1, efc, emc, tbr, lseg, sp3;
   int va, v3, ok1, nh;
   w.SetBranchAddress("hs", &hs); w.SetBranchAddress("bA", &bA); w.SetBranchAddress("br", &br); w.SetBranchAddress("b1", &b1);
   w.SetBranchAddress("efc", &efc); w.SetBranchAddress("emc", &emc); w.SetBranchAddress("tbr", &tbr);
   w.SetBranchAddress("lseg", &lseg); w.SetBranchAddress("sp3", &sp3); w.SetBranchAddress("va", &va);
   w.SetBranchAddress("v3", &v3); w.SetBranchAddress("ok1", &ok1); w.SetBranchAddress("nh", &nh);
   const int nm = 4, nb = 7;
   const double bl[nb + 1] = {0.20, 0.25, 0.30, 0.35, 0.45, 0.60, 0.75, 0.90};
   const char *name[nm] = {"Spyral", "3D", "3D hit-refit + arc angle (adopted)", "species-blind EOM"};
   const int col[nm] = {kRed + 1, kGray + 2, kBlue + 2, kGreen + 2};
   std::vector<double> bias[nm][nb], land[nm][nb], mass[nb];
   std::vector<long> nTrue(nb, 0);
   for (Long64_t i = 0; i < w.GetEntries(); ++i) {
      w.GetEntry(i);
      if (!va || !v3 || !ok1 || sp3 >= 6 || nh < 30 || tbr <= 0) continue;
      int st = -1;
      for (int k = 0; k < nb; ++k) if (tbr > bl[k] && tbr < bl[k + 1]) st = k;
      double B[nm] = {br, b1, bA, efc};
      for (int m = 0; m < nm; ++m) {
         if (B[m] <= 0) continue;
         if (st >= 0) bias[m][st].push_back(std::log(B[m] / tbr));
         int sr = -1;
         for (int k = 0; k < nb; ++k) if (B[m] > bl[k] && B[m] < bl[k + 1]) sr = k;
         if (sr >= 0 && hs > 0 && lseg > 5)
            land[m][sr].push_back(std::log(hs / std::sqrt(eom::meanS(tab[0], B[m] * 299.792458, eom::massMeV(1), lseg) / ksim)));
      }
      if (st >= 0) { ++nTrue[st]; if (emc > 0) mass[st].push_back(std::log(emc / eom::massMeV(2))); }
   }
   auto stat = [](std::vector<double> &v, double &med, double &hw) {
      if (v.size() < 10) return false;
      std::sort(v.begin(), v.end());
      size_t n = v.size();
      med = v[n / 2]; hw = 0.5 * (v[n * 84 / 100] - v[n * 16 / 100]);
      return true;
   };
   auto *c = new TCanvas("c", "c", 1800, 1200);
   c->Divide(2, 2);
   const char *ttl[4] = {"(a) B#rho bias, tritons (sim);true B#rho [T m];median ln(B#rho / B#rho_{true})",
                         "(b) B#rho resolution;true B#rho [T m];half 68 % width of ln(B#rho / B#rho_{true})",
                         "(c) where tritons land: ln(#sqrt{dE/dx} / deuteron prediction) at the method's own B#rho;reconstructed B#rho [T m];#delta from the deuteron locus",
                         "(d) species-blind fitted mass, tritons;true B#rho [T m];ln(m_{fit} / m_{t})"};
   const double yr[4][2] = {{-0.5, 0.1}, {0, 0.3}, {-0.1, 0.8}, {-0.6, 0.4}};
   printf("tritons per true-Brho slice:");
   for (int k = 0; k < nb; ++k) printf(" %.2f-%.2f:%ld", bl[k], bl[k + 1], nTrue[k]);
   printf("\n");
   for (int p = 0; p < 4; ++p) {
      c->cd(p + 1); gPad->SetGrid();
      auto *fr = gPad->DrawFrame(0.18, yr[p][0], 0.92, yr[p][1]);
      fr->SetTitle(ttl[p]);
      if (p == 2) {
         auto *band = new TBox(0.18, -0.08, 0.92, 0.08); band->SetFillColorAlpha(kOrange, 0.35); band->Draw();
         auto *t = new TLatex(0.2, 0.09, "deuterons (data half-width)"); t->SetTextSize(0.035); t->Draw();
      }
      auto *lg = p == 2 ? new TLegend(0.35, 0.70, 0.89, 0.89) : new TLegend(0.35, 0.13, 0.89, 0.33);
      for (int m = 0; m < (p == 3 ? 1 : nm); ++m) {
         auto *g = new TGraphErrors();
         for (int k = 0; k < nb; ++k) {
            double med, hw;
            auto &v = p == 3 ? mass[k] : (p == 2 ? land[m][k] : bias[m][k]);
            if (!stat(v, med, hw)) continue;
            double x = 0.5 * (bl[k] + bl[k + 1]) + (m - 1.5) * 0.006;
            int n = g->GetN();
            if (p == 1) g->SetPoint(n, x, hw); else { g->SetPoint(n, x, med); g->SetPointError(n, 0, hw); }
         }
         g->SetMarkerStyle(20 + m); g->SetMarkerSize(1.4);
         int cc = p == 3 ? col[3] : col[m];
         g->SetMarkerColor(cc); g->SetLineColor(cc); g->SetLineWidth(2);
         g->Draw(p == 1 ? "LP" : "P");
         lg->AddEntry(g, p == 3 ? name[3] : name[m], "lp");
      }
      if (p == 3) { auto *l = new TLine(0.18, -0.405, 0.92, -0.405); l->SetLineStyle(2); l->Draw(); auto *t = new TLatex(0.2, -0.38, "deuteron mass"); t->SetTextSize(0.035); t->Draw(); }
      if (p == 0) { auto *l = new TLine(0.18, 0, 0.92, 0); l->SetLineStyle(2); l->Draw(); }
      lg->Draw();
   }
   // table
   printf("%-36s", "true Brho slice");
   for (int k = 0; k < nb; ++k) printf("   %.2f-%.2f  ", bl[k], bl[k + 1]);
   printf("\n");
   for (int m = 0; m < nm; ++m) {
      printf("bias  %-30s", name[m]);
      for (int k = 0; k < nb; ++k) { double a, b; if (stat(bias[m][k], a, b)) printf(" %+.3f[%.3f]", a, b); else printf("      --     "); }
      printf("\n");
   }
   for (int m = 0; m < nm; ++m) {
      printf("land  %-30s", name[m]);
      for (int k = 0; k < nb; ++k) { double a, b; if (stat(land[m][k], a, b)) printf(" %+.3f[%.3f]", a, b); else printf("      --     "); }
      printf("\n");
   }
   printf("mass  %-30s", "species-blind");
   for (int k = 0; k < nb; ++k) { double a, b; if (stat(mass[k], a, b)) printf(" %+.3f[%.3f]", a, b); else printf("      --     "); }
   printf("\n");
   c->SaveAs(png);
}
