/// @file fwd_dwba_C14.C
/// @brief Compare each level's FORWARD angular distribution (default 20-60 deg) with the DWBA.
///
/// The full 25-135 deg fit asks the collective model to describe angles where it is least reliable:
/// backward of ~70 deg the cross section is small and whatever compound or multi-step strength there
/// is contributes a larger FRACTION of it. Forward of 60 deg the reaction is most nearly direct, and
/// it is also where L = 2 and L = 3 separate by slope (over 30->60 deg the normalised Perey shapes
/// fall to 0.75 and 0.61 respectively) rather than by the backward rise.
///
/// Reports, per level, chi2/ndf against each L, and for the 7.012 the B(E2) implied by the forward
/// region alone -- to be compared with the same number from the full range.
///
///   root -b -q 'fwd_dwba_C14.C("plots/fit_angles_ps_dist_tcfwd.root",20,60)'
namespace fwd {
TGraph *rd(TString f)
{
   auto *g = new TGraph();
   std::ifstream in(f.Data());
   double a, b;
   while (in >> a >> b) if (b > 0) g->SetPoint(g->GetN(), a, b);
   return g;
}
/// The .beta file is "LX betaNuclear betaCoulomb B(EL)[e2 barn^LX] Rdeformation" -- B(EL) is the
/// FOURTH column. Reading the first returns LX (= 2), which is 20x the real BELX of 0.1.
double belxOf(TString f)
{
   std::ifstream bi(f.Data()); std::string h; std::getline(bi, h);
   int lx = 0; double bN = 0, bC = 0, belx = 0, Rd = 0;
   bi >> lx >> bN >> bC >> belx >> Rd;
   return belx;
}
} // namespace fwd

void fwd_dwba_C14(TString distFile = "plots/fit_angles_ps_dist_tcfwd.root",
                  Double_t loA = 20, Double_t hiA = 60, TString pot = "P")
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TString pdir = here + "/../ptolemy/dat/";
   TFile *f = TFile::Open(here + "/" + distFile);
   if (!f || f->IsZombie()) { printf("\033[1;31mno %s\033[0m\n", distFile.Data()); return; }

   const char *lv[5]   = {"lvl0", "lvl1", "lvl2", "lvl3", "lvl4"};
   const char *nm[5]   = {"6.091 1-", "6.728 3-", "7.012 2+", "7.341 2-", "8.317 2+"};
   const char *ltag[5] = {"6094_1m", "6728_3m", "7012_2p", "7012_2p", "8317_2p"};

   printf("\n  forward-region test, %.0f-%.0f deg, potential %s\n", loA, hiA, pot.Data());
   printf("  chi2/ndf of a single free scale against each multipole\n\n");
   printf("  %-11s %5s |%s | %s\n", "level", "npts",
          "     L=0     L=1     L=2     L=3", "best");
   printf("  %s\n", TString('-', 62).Data());

   std::vector<TGraph *> B(4);
   for (int L = 0; L < 4; ++L) B[L] = fwd::rd(Form("%smdaAll_L%d.dat", pdir.Data(), L));

   auto *cv = new TCanvas("cfwd", "", 1500, 900); cv->Divide(3, 2);
   for (int i = 0; i < 5; ++i) {
      auto *g = (TGraphErrors *)f->Get(lv[i]);
      if (!g) continue;
      std::vector<double> th, y, e;
      for (int p = 0; p < g->GetN(); ++p) {
         double x, v; g->GetPoint(p, x, v); double er = g->GetErrorY(p);
         if (x < loA || x > hiA || er <= 0 || v <= 0) continue;
         th.push_back(x); y.push_back(v); e.push_back(er);
      }
      if (th.size() < 3) { printf("  %-11s  too few points in range\n", nm[i]); continue; }
      double c2[4], amp[4];
      for (int L = 0; L < 4; ++L) {
         double num = 0, den = 0;
         for (size_t k = 0; k < th.size(); ++k) {
            double fv = B[L]->Eval(th[k]);
            num += y[k] * fv / (e[k] * e[k]); den += fv * fv / (e[k] * e[k]);
         }
         amp[L] = std::max(0.0, num / den);
         double s = 0;
         for (size_t k = 0; k < th.size(); ++k) s += std::pow((y[k] - amp[L] * B[L]->Eval(th[k])) / e[k], 2);
         c2[L] = s / (th.size() - 1);
      }
      int best = 0; for (int L = 1; L < 4; ++L) if (c2[L] < c2[best]) best = L;
      printf("  %-11s %5zu |", nm[i], th.size());
      for (int L = 0; L < 4; ++L) printf(" %7.2f", c2[L]);
      printf(" | L=%d\n", best);

      cv->cd(i + 1); gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
      auto *gd = new TGraphErrors((int)th.size());
      for (size_t k = 0; k < th.size(); ++k) { gd->SetPoint(k, th[k], y[k]); gd->SetPointError(k, 0, e[k]); }
      gd->SetMarkerStyle(20); gd->SetMarkerSize(1.2);
      gd->SetTitle(Form("%s   (%.0f-%.0f deg);#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]", nm[i], loA, hiA));
      gd->Draw("AP");
      auto *lg = new TLegend(0.55, 0.66, 0.93, 0.90); lg->SetBorderSize(0); lg->SetFillStyle(0);
      lg->SetTextSize(0.038);
      const int cl[4] = {kGray + 2, kOrange + 7, kAzure + 2, kGreen + 3};
      for (int L = 0; L < 4; ++L) {
         if (L == 0 || L == 1) continue;      // draw only the two that compete here
         auto *q = new TGraph((int)th.size());
         for (size_t k = 0; k < th.size(); ++k) q->SetPoint(k, th[k], amp[L] * B[L]->Eval(th[k]));
         q->SetLineColor(cl[L]); q->SetLineWidth(3); q->SetLineStyle(L == best ? 1 : 2); q->Draw("L same");
         lg->AddEntry(q, Form("L = %d   #chi^{2}/ndf %.2f", L, c2[L]), "l");
      }
      lg->Draw();
   }

   // ---- B(E2) of the 7.012 from the forward region alone -------------------------------------
   auto *g2 = (TGraphErrors *)f->Get("lvl2");
   auto *gp = fwd::rd(pdir + "omp_" + pot + "_7012_2p.dat");
   double belx = fwd::belxOf(pdir + "omp_" + pot + "_7012_2p.beta");
   if (g2 && gp->GetN() && belx > 0) {
      double num = 0, den = 0; int n = 0;
      for (int p = 0; p < g2->GetN(); ++p) {
         double x, v; g2->GetPoint(p, x, v); double er = g2->GetErrorY(p);
         if (x < loA || x > hiA || er <= 0 || v <= 0) continue;
         double fv = gp->Eval(x);
         num += v * fv / (er * er); den += fv * fv / (er * er); ++n;
      }
      double k = den > 0 ? num / den : 0;
      const double WU = 0.05940 * std::pow(14.0, 4.0 / 3.0);
      printf("\n  7.012 B(E2) from %0.f-%0.f deg alone, %s curve (BELX %.3g e2b^2), %d points:\n",
             loA, hiA, pot.Data(), belx, n);
      printf("    B(E2)up = %.1f e2fm4   = %.2f W.u. up   = %.2f W.u. down\n",
             k * belx * 1e4, k * belx * 1e4 / WU, k * belx * 1e4 / WU / 5);
      printf("    (full 25-135 deg gave 32.0 e2fm4 = 3.19 W.u. down; ENSDF 1.80)\n");
   }
   cv->cd(6);
   auto *pt = new TPaveText(0.03, 0.1, 0.97, 0.9); pt->SetBorderSize(0); pt->SetFillStyle(0);
   pt->SetTextAlign(12);
   pt->AddText(Form("#bf{forward region %.0f-%.0f deg}", loA, hiA));
   pt->AddText("");
   pt->AddText("dashed = the losing multipole");
   pt->AddText("solid  = the preferred one");
   pt->Draw();
   TString out = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";
   cv->SaveAs(out + "24_forward_dwba.png");
   printf("\n  wrote %s24_forward_dwba.png\n\n", out.Data());
}
