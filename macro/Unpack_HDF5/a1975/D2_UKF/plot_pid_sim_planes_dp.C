/// @file plot_pid_sim_planes_dp.C
/// @brief PID planes of the simulated (d,t) tritons (g.s. + ex2) for each Brho method (2026-11-01): Brho on y, the
///        3D window sqrt(dE/dx) hs on x (same for every panel, so only the Brho differs). Curves: CATIMA prediction
///        sqrt(meanS(p, m, lseg)/k_sim) for deuterons (no deuterons in this sim) and tritons at lseg = 170 mm (typical),
///        drawn down to the Brho where that ion stops within 170 mm.
///   root -l -b -q 'plot_pid_sim_planes_dp.C("plots/pid_sim_planes.png")'
#include "eom_brho_dp.h"
void plot_pid_sim_planes_dp(TString png, double ksim = 1.982e-5, double lsegRef = 170,
                            TString files = "/mnt/f/a1975/pid3d_sbpid/sim_gs.root,/mnt/f/a1975/pid3d_sbpid/sim_ex2.root")
{
   gSystem->Load("libAtTools.so");
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kViridis);
   eom::Table tab[3];
   eom::makeTables(tab);
   TChain w("w");
   for (auto *o : *files.Tokenize(",")) w.Add(((TObjString *)o)->GetString());
   TString sel = "va&&v3&&ok1&&sp3<6&&nh>=30&&hs>0&&tbr>0";
   const char *B[5] = {"tbr", "br", "b1", "bA", "efc"};
   const char *nm[5] = {"truth B#rho", "Spyral B#rho", "3D B#rho", "3D hit-refit + arc (adopted)", "species-blind EOM B#rho"};
   w.Draw("hs", sel, "goff");
   std::vector<double> x(w.GetV1(), w.GetV1() + w.GetSelectedRows());
   std::sort(x.begin(), x.end());
   const double xhi = x[x.size() * 995 / 1000] * 1.3;
   // prediction curves
   TGraph *gp[2];
   for (int s = 0; s < 2; ++s) {
      gp[s] = new TGraph();
      // from high Brho down; stop where the ion would stop inside lsegRef (the mean then falls, curve turns back)
      double xprev = 0;
      for (double b = 1.0; b >= 0.12; b -= 0.005) {
         double ms = eom::meanS(tab[0], b * 299.792458, eom::massMeV(s + 1), lsegRef);
         if (ms <= 0 || std::sqrt(ms / ksim) < xprev) break;
         xprev = std::sqrt(ms / ksim);
         gp[s]->SetPoint(gp[s]->GetN(), xprev, b);
      }
      gp[s]->SetLineWidth(2); gp[s]->SetLineStyle(2); gp[s]->SetLineColor(s == 0 ? kOrange + 7 : kRed + 1);
   }
   auto *c = new TCanvas("c", "c", 2000, 450);
   c->Divide(5, 1);
   for (int m = 0; m < 5; ++m) {
      c->cd(m + 1);
      gPad->SetLogz(); gPad->SetGrid(); gPad->SetRightMargin(0.12);
      auto *h = new TH2D(Form("h%d", m), Form("sim tritons: %s;#sqrt{dE/dx} (3D window, hs);B#rho [T m]", nm[m]), 200, 0, xhi, 200, 0, 1.0);
      w.Draw(Form("%s:hs>>h%d", B[m], m), sel + Form("&&%s>0", B[m]), "goff");
      h->Draw("colz");
      gp[0]->Draw("L"); gp[1]->Draw("L");
      if (m == 0) {
         auto *lg = new TLegend(0.40, 0.75, 0.87, 0.88);
         lg->AddEntry(gp[0], "deuteron prediction", "l"); lg->AddEntry(gp[1], "triton prediction", "l"); lg->Draw();
      }
   }
   c->SaveAs(png);
}
