/// @file plot_pid_gf_sim_dp.C
/// @brief PID from genfit Brho + predicted energy loss, on the d + t sims (2026-11-01). Each track is fitted by genfit as
///        d and as t (gf_hyp_run_dp.sh, measSigma 0.6). For each hypothesis h the vertex Brho B_h predicts the mean
///        window dE/dx (CATIMA, eom::meanS over the track's own PID segment lseg, k_sim); residual
///        r_h = ln(hs_measured / hs_predicted_h). Panels: (1) hs vs genfit d-hypothesis Brho, (2) hs vs t-hypothesis
///        Brho, both with the CATIMA d and t curves (lseg 170 mm), (3) r_t vs r_d: a deuteron sits at r_d ~ 0, a
///        triton at r_t ~ 0; the dashed diagonals |r_d| = |r_t| are the decision boundary. Sim d blue, sim t red.
///   root -l -b -q 'plot_pid_gf_sim_dp.C("plots/pid_gf_sim.png")'
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/eom_brho_dp.h"
/// r_d, r_t per selected track (kept out of the main function: inlined there, cling miscompiled the CATIMA table
/// build at the top of the function -- every table entry came back NaN, 2026-11-01)
void residuals(TTree *T, TString sel, eom::Table *tab, double k, std::vector<double> &X, std::vector<double> &Y)
{
   Long64_t N = T->Draw("hs:g.bx_d:g.bx_t:lseg", sel, "goff");
   std::vector<double> hs(T->GetV1(), T->GetV1() + N), bd(T->GetV2(), T->GetV2() + N), bt(T->GetV3(), T->GetV3() + N),
      ls(T->GetV4(), T->GetV4() + N);
   for (Long64_t i = 0; i < N; ++i) {
      double pd = eom::meanS(tab[0], bd[i] * 299.792458, eom::massMeV(1), ls[i]);
      double pt = eom::meanS(tab[0], bt[i] * 299.792458, eom::massMeV(2), ls[i]);
      if (!(pd > 0) || !(pt > 0)) continue;
      X.push_back(std::log(hs[i] / std::sqrt(pd / k))); Y.push_back(std::log(hs[i] / std::sqrt(pt / k)));
   }
}
void plot_pid_gf_sim_dp(TString png, TString G = "/mnt/f/a1975/gf_hyp/s06/", double k = 1.982e-5)
{
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   eom::Table tab[3];
   eom::makeTables(tab);
   TChain d("w"), t("w"), dg("g"), tg("g");
   for (int s = 5001; s <= 5004; ++s) { d.Add(Form("/mnt/f/a1975_C16_dd_sim/pv_el_s%d.root", s)); dg.Add(Form("%shyp_el_s%d.root", G.Data(), s)); }
   for (const char *s : {"gs", "ex2"}) { t.Add(Form("/mnt/f/a1975/pid3d_sbpid/sim_%s.root", s)); tg.Add(Form("%shyp_%s.root", G.Data(), s)); }
   d.AddFriend(&dg, "g"); t.AddFriend(&tg, "g");
   TTree *T[2] = {&d, &t};
   TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30&&hs>0&&lseg>5&&g.ok_d&&g.ok_t";
   auto *c = new TCanvas("c", "c", 1950, 650);
   c->Divide(3, 1);
   const int col[2] = {kBlue + 1, kRed + 1};
   TGraph *cur[2];
   for (int s = 0; s < 2; ++s) {
      cur[s] = new TGraph();
      double xprev = 0;
      for (double b = 1.0; b >= 0.12; b -= 0.005) {
         double ms = eom::meanS(tab[0], b * 299.792458, eom::massMeV(s + 1), 170);
         if (ms <= 0 || std::sqrt(ms / k) < xprev) break;
         xprev = std::sqrt(ms / k);
         cur[s]->SetPoint(cur[s]->GetN(), xprev, b);
      }
      cur[s]->SetLineWidth(2); cur[s]->SetLineStyle(2); cur[s]->SetLineColor(s ? kRed + 2 : kBlue + 3);
   }
   for (int p = 0; p < 2; ++p) {
      c->cd(p + 1); gPad->SetGrid();
      auto *fr = gPad->DrawFrame(0, 0, 26, 1.0);
      fr->SetTitle(Form("genfit B#rho, fitted as %s;#sqrt{dE/dx} (3D window, hs);B#rho_{genfit} [T m]", p ? "TRITON" : "DEUTERON"));
      for (int s = 0; s < 2; ++s) {
         long n = T[s]->Draw(Form("g.bx_%s:hs", p ? "t" : "d"), sel, "goff");
         auto *g = new TGraph(n, T[s]->GetV2(), T[s]->GetV1());
         g->SetMarkerStyle(1); g->SetMarkerColorAlpha(col[s], 0.6); g->Draw("P");
      }
      cur[0]->Draw("L"); cur[1]->Draw("L");
      auto *lg = new TLegend(0.45, 0.72, 0.88, 0.88);
      auto *gd = new TGraph(); gd->SetMarkerStyle(20); gd->SetMarkerColor(col[0]);
      auto *gt = new TGraph(); gt->SetMarkerStyle(20); gt->SetMarkerColor(col[1]);
      lg->AddEntry(gd, "sim deuterons", "p"); lg->AddEntry(gt, "sim tritons", "p");
      lg->AddEntry(cur[0], "CATIMA d (170 mm)", "l"); lg->AddEntry(cur[1], "CATIMA t (170 mm)", "l"); lg->Draw();
   }
   c->cd(3); gPad->SetGrid();
   auto *fr = gPad->DrawFrame(-0.6, -0.6, 0.6, 0.6);
   fr->SetTitle("residual to each hypothesis, B#rho_{genfit} < 0.6;r_{d} = ln(hs / hs_{pred, d});r_{t} = ln(hs / hs_{pred, t})");
   for (int s = 0; s < 2; ++s) {
      std::vector<double> X, Y;
      residuals(T[s], sel + "&&g.bx_d<0.6", tab, k, X, Y);
      auto *g = new TGraph(X.size(), X.data(), Y.data());
      g->SetMarkerStyle(7); g->SetMarkerColorAlpha(col[s], 0.5); g->Draw("P");
   }
   for (int sg : {1, -1}) { auto *l = new TLine(-0.6, -0.6 * sg, 0.6, 0.6 * sg); l->SetLineStyle(2); l->Draw(); }
   c->SaveAs(png);
}
