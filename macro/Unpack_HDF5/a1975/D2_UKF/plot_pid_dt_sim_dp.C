/// @file plot_pid_dt_sim_dp.C
/// @brief Simulated deuterons (elastic d(16C,d), C16_dd_el_sim.C) AND tritons ((d,t) g.s. + ex2) through each Brho
///        method, same gain, same 3D window dE/dx hs (2026-11-01). Top: PID planes (Brho on y), d blue, t red.
///        Bottom: per slice of each method's own Brho, the d/t separation S = (median_t - median_d) / (hw_t + hw_d) of
///        ln(hs * Brho) (hw = half 68 % width), and the misidentified fractions at the midpoint.
///        Trees: pid_sbpid_dp.C on the sims + spiral_fit_dp.C (same rows, friend "s").
///   root -l -b -q 'plot_pid_dt_sim_dp.C("plots/pid_dt_sim.png")'
void plot_pid_dt_sim_dp(TString png, TString spB = "s.b21")
{
   gStyle->SetOptStat(0);
   TChain d("w"), t("w"), ds("w"), tsp("w");
   for (int s = 5001; s <= 5004; ++s) {
      d.Add(Form("/mnt/f/a1975_C16_dd_sim/pv_el_s%d.root", s)); ds.Add(Form("/mnt/f/a1975_C16_dd_sim/sp_el_s%d.root", s));
   }
   for (const char *s : {"gs", "ex2"}) {
      t.Add(Form("/mnt/f/a1975/pid3d_sbpid/sim_%s.root", s)); tsp.Add(Form("/mnt/f/a1975/pid3d_spiral/sim_%s.root", s));
   }
   d.AddFriend(&ds, "s"); t.AddFriend(&tsp, "s");
   TTree *T[2] = {&d, &t};
   const int nm = 4;
   TString B[nm] = {"tbr", "bA", "br", spB};
   const char *name[nm] = {"truth B#rho", "circle + arc (adopted)", "Spyral B#rho", "spiral (shrinking radius)"};
   TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30&&hs>0&&tbr>0";
   auto *c = new TCanvas("c", "c", 2000, 1150);
   auto *top = new TPad("top", "", 0, 0.42, 1, 1); top->Draw();
   auto *bot = new TPad("bot", "", 0, 0, 1, 0.42); bot->Draw();
   top->Divide(nm, 1);
   for (int m = 0; m < nm; ++m) {
      top->cd(m + 1); gPad->SetGrid();
      auto *fr = gPad->DrawFrame(0, 0, 26, 1.0);
      fr->SetTitle(Form("%s;#sqrt{dE/dx} (3D window, hs);B#rho [T m]", name[m]));
      for (int sp = 0; sp < 2; ++sp) {
         long n = T[sp]->Draw(Form("%s:hs", B[m].Data()), sel + Form("&&%s>0", B[m].Data()), "goff");
         auto *g = new TGraph(n, T[sp]->GetV2(), T[sp]->GetV1());
         g->SetMarkerStyle(1); g->SetMarkerColorAlpha(sp == 0 ? kBlue + 1 : kRed + 1, 0.6);
         g->Draw("P");
      }
      if (m == 0) {
         auto *lg = new TLegend(0.55, 0.78, 0.88, 0.88);
         auto *gd = new TGraph(); gd->SetMarkerStyle(20); gd->SetMarkerColor(kBlue + 1);
         auto *gt = new TGraph(); gt->SetMarkerStyle(20); gt->SetMarkerColor(kRed + 1);
         lg->AddEntry(gd, "deuterons (sim)", "p"); lg->AddEntry(gt, "tritons (sim)", "p"); lg->Draw();
      }
   }
   // separation vs Brho
   const double bl[] = {0.20, 0.25, 0.30, 0.35, 0.45, 0.60, 0.75, 0.90};
   const int nb = 7;
   bot->Divide(2, 1);
   TGraph *gS[nm], *gM[nm];
   const int col[nm] = {kBlack, kBlue + 2, kRed + 1, kGreen + 2};
   printf("%-28s", "own-Brho slice");
   for (int k = 0; k < nb; ++k) printf("    %.2f-%.2f     ", bl[k], bl[k + 1]);
   printf("\n");
   for (int m = 0; m < nm; ++m) {
      gS[m] = new TGraph(); gM[m] = new TGraph();
      printf("%-28s", name[m]);
      for (int k = 0; k < nb; ++k) {
         std::vector<double> v[2];
         for (int sp = 0; sp < 2; ++sp) {
            long n = T[sp]->Draw(Form("log(hs*%s)", B[m].Data()),
                                 sel + Form("&&%s>%g&&%s<%g", B[m].Data(), bl[k], B[m].Data(), bl[k + 1]), "goff");
            v[sp].assign(T[sp]->GetV1(), T[sp]->GetV1() + n);
            std::sort(v[sp].begin(), v[sp].end());
         }
         if (v[0].size() < 15 || v[1].size() < 15) { printf("        --         "); continue; }
         auto q = [](std::vector<double> &x, double f) { return x[std::min(x.size() - 1, (size_t)(f * x.size()))]; };
         double md = q(v[0], 0.5), mt = q(v[1], 0.5);
         double hd = 0.5 * (q(v[0], 0.84) - q(v[0], 0.16)), ht = 0.5 * (q(v[1], 0.84) - q(v[1], 0.16));
         double S = (mt - md) / (hd + ht);
         // misidentification at the point dividing the gap in proportion to the widths
         double cut = md + (mt - md) * hd / (hd + ht);
         double fd = double(v[0].end() - std::upper_bound(v[0].begin(), v[0].end(), cut)) / v[0].size();
         double ft = double(std::upper_bound(v[1].begin(), v[1].end(), cut) - v[1].begin()) / v[1].size();
         double x = 0.5 * (bl[k] + bl[k + 1]);
         gS[m]->SetPoint(gS[m]->GetN(), x, S);
         gM[m]->SetPoint(gM[m]->GetN(), x, 100 * ft);
         printf("  S %4.2f t>d %4.1f%% d>t %4.2f%%", S, 100 * ft, 100 * fd);
      }
      printf("\n");
   }
   for (int p = 0; p < 2; ++p) {
      bot->cd(p + 1); gPad->SetGrid();
      auto *fr = gPad->DrawFrame(0.18, 0, 0.92, p == 0 ? 4.0 : 60);
      fr->SetTitle(p == 0 ? "d/t separation (#Delta median / sum of half-widths), sim;method's own B#rho [T m];separation"
                          : "tritons that fall on the deuteron side of the cut, sim;method's own B#rho [T m];t #rightarrow d [%]");
      auto *lg = new TLegend(0.5, 0.65, 0.89, 0.89);
      for (int m = 0; m < nm; ++m) {
         TGraph *g = p == 0 ? gS[m] : gM[m];
         g->SetMarkerStyle(20 + m); g->SetMarkerSize(1.4); g->SetMarkerColor(col[m]); g->SetLineColor(col[m]); g->SetLineWidth(2);
         g->Draw("LP");
         lg->AddEntry(g, name[m], "lp");
      }
      lg->Draw();
   }
   c->SaveAs(png);
}
