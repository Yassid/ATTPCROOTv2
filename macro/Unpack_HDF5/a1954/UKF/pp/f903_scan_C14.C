/// @file f903_scan_C14.C
/// @brief How much can the 6.903 (0-) inflate the "7.012" yield, and hence B(E2)?
///
/// The 6.903 and 7.012 sit 0.109 MeV apart against a resolution sigma of 0.143 -- 0.76 sigma,
/// which is ONE peak. Fitted as a free component the model degenerates (6.728 collapses to zero
/// beyond 65 deg and 6.903 takes its strength). So the 6.903 amplitude is instead TIED to a fixed
/// fraction f of the 7.012 one -- no free parameter -- and f is scanned. The result is a band on
/// B(E2), not a point, which is the honest form given the resolution.
///
///   root -b -q 'f903_scan_C14.C()'
void f903_scan_C14(TString curve = "exc_7012_2p", TString tagPrefix = "f",
                   Double_t fitLo = 25, Double_t fitHi = 135)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TString pdir = here + "/../ptolemy/dat/";
   const int NF = 5;
   const double F[NF] = {0.0, 0.10, 0.20, 0.30, 0.50};
   const char *SUF[NF] = {"00", "010", "020", "030", "050"};
   std::vector<TString> TGv;
   for (int i = 0; i < NF; ++i) TGv.push_back(tagPrefix + SUF[i]);

   // the Ptolemy 7.012 curve and the BELX it was run at
   auto *gp = new TGraph();
   { std::ifstream in((pdir + curve + ".dat").Data()); double a, b;
     while (in >> a >> b) if (b > 0) gp->SetPoint(gp->GetN(), a, b); }
   double belx = 0, bN, bC, Rd; int lx;
   { std::ifstream bi((pdir + curve + ".beta").Data()); std::string h; std::getline(bi, h);
     bi >> lx >> bN >> bC >> belx >> Rd; }
   if (!gp->GetN() || belx <= 0) { printf("\033[1;31mno Ptolemy curve %s\033[0m\n", curve.Data()); return; }
   printf("\n  DWBA curve: %s   (BELX %.3g e2 b^L)\n", curve.Data(), belx);

   const double WU = 0.05940 * std::pow(14.0, 4.0 / 3.0);   // e2 fm4
   const double LIT = 5 * 1.8 * WU;                          // ENSDF 1.8(3) W.u. DOWN -> UP
   printf("\n  literature B(E2)up = %.1f e2fm4 (ENSDF 1.8(3) W.u. downward, x5)\n", LIT);
   printf("\n  %6s %10s %12s %10s %10s %8s %10s\n",
          "f903", "sum sigma", "B(E2)up", "W.u.", "ratio EM", "rms", "d_h/d_EM");
   std::vector<double> fv, bv, rv, sv;
   for (int k = 0; k < NF; ++k) {
      TFile *f = TFile::Open(here + "/plots/fit_angles_ps_dist_" + TGv[k] + ".root");
      auto *gd = f && !f->IsZombie() ? (TGraphErrors *)f->Get("lvl2") : nullptr;   // 7.012
      if (!gd) { printf("  %6.2f   missing %s\n", F[k], TGv[k].Data()); continue; }
      double s = 0, tot = 0; int n = 0;
      for (int j = 0; j < gd->GetN(); ++j) {
         double x = gd->GetX()[j], y = gd->GetY()[j];
         if (x < fitLo || x > fitHi || y <= 0) continue;
         double p = gp->Eval(x);
         if (p > 0) { s += std::log(y / p); ++n; tot += y; }
      }
      if (!n) continue;
      double kk = std::exp(s / n), r = 0; int m = 0;
      for (int j = 0; j < gd->GetN(); ++j) {
         double x = gd->GetX()[j], y = gd->GetY()[j];
         if (x < fitLo || x > fitHi || y <= 0) continue;
         double p = kk * gp->Eval(x);
         if (p > 0) { r += std::pow(std::log(y / p), 2); ++m; }
      }
      double B = kk * belx * 1e4;              // e2 b2 -> e2 fm4
      double rms = std::sqrt(r / m), rat = B / LIT;
      printf("  %6.2f %10.2f %12.1f %10.1f %10.2f %8.3f %10.2f\n",
             F[k], tot, B, B / WU, rat, rms, std::sqrt(rat));
      fv.push_back(F[k]); bv.push_back(B); rv.push_back(rms); sv.push_back(tot);
   }
   if (fv.size() < 2) return;
   printf("\n  B(E2) falls %.0f%% from f903 = 0 to %.2f\n", 100 * (1 - bv.back() / bv.front()), fv.back());
   printf("  the shape quality (rms) changes by %.3f over that range -- the data do NOT\n"
          "  prefer one f903, which is exactly why this is a band and not a measurement\n", 
          *std::max_element(rv.begin(), rv.end()) - *std::min_element(rv.begin(), rv.end()));

   auto *c = new TCanvas("cf9", "", 1200, 500); c->Divide(2, 1);
   c->cd(1); gPad->SetGridx(); gPad->SetGridy();
   auto *g1 = new TGraph(fv.size(), &fv[0], &bv[0]);
   // ROOT splits SetTitle on ";" -- the one inside "B(E2;0+->2+)" was being read as the axis
   // separator, so the plot came out titled "B(E2" with every axis label shifted by one: the
   // y axis carried the x label and the real y label was dropped. No semicolon in the title text.
   g1->SetTitle("B(E2) 0^{+}#rightarrow2^{+} vs the tied 6.903 fraction;"
                "f_{903} = A(6.903)/A(7.012);B(E2)#uparrow [e^{2}fm^{4}]");
   g1->SetMarkerStyle(20); g1->SetMarkerSize(1.6); g1->SetLineWidth(3); g1->SetMinimum(0);
   g1->SetMaximum(1.2 * (*std::max_element(bv.begin(), bv.end()))); g1->Draw("ALP");
   auto *lb = new TLine(fv.front(), LIT, fv.back(), LIT);
   lb->SetLineColor(kRed + 1); lb->SetLineWidth(3); lb->Draw();
   auto *bx = new TBox(fv.front(), LIT - 5 * 0.3 * WU, fv.back(), LIT + 5 * 0.3 * WU);
   bx->SetFillColorAlpha(kRed + 1, 0.20); bx->Draw();
   TLatex t; t.SetNDC(); t.SetTextSize(0.042); t.SetTextColor(kRed + 1);
   t.DrawLatex(0.45, 0.26, "ENSDF, 1.8(3) W.u. #times 5");
   c->cd(2); gPad->SetGridx(); gPad->SetGridy();
   auto *g2 = new TGraph(fv.size(), &fv[0], &rv[0]);
   g2->SetTitle("shape quality -- flat means the data cannot choose f_{903};f_{903};rms ln(data/Ptolemy)");
   g2->SetMarkerStyle(21); g2->SetMarkerSize(1.6); g2->SetLineWidth(3);
   g2->SetMinimum(0); g2->SetMaximum(1.0); g2->Draw("ALP");
   TString out = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";
   c->SaveAs(out + "10_f903_scan_" + tagPrefix + ".png");
   printf("\n  wrote %s10_f903_scan_%s.png\n\n", out.Data(), tagPrefix.Data());
}
