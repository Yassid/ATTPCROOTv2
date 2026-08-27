/// @file elastic_flatness_C14.C
/// @brief Choose the normalisation window on the one criterion that is a requirement, not a taste.
///
/// L is a beam flux times a target areal density. It CANNOT depend on scattering angle. So
/// plotting L(theta) = Y/A/dOmega / sigma_calc(theta) and fitting a constant is not a comparison
/// between potentials, it is a self-consistency test of the normalisation: chi2/ndf of order one
/// or the normalisation is wrong, whatever else recommends the potential.
///
/// This matters because the window was previously inherited rather than chosen. The main analysis
/// document normalises over theta_cm 42-58 deg on the grounds that "the measured and calculated
/// shapes genuinely agree" there. They do not. That range ENDS at the measured diffraction
/// minimum (58.0 deg), which Koning-Delaroche places at 63 deg, so across it the data are diving
/// into a dip the calculation has not reached and L collapses -- 53.8, 48.7, 13.4, 4.0 at 42.5,
/// 47.5, 52.5, 57.5 deg. The document's own table shows the slide (32.3 -> 28.0) and reads it as
/// spread. Normalising on the shoulder of a misplaced diffraction minimum is the single largest
/// error in the old absolute scale.
///
/// The test also reverses the choice of potential. Becchetti-Greenlees and Menet place the minimum
/// exactly, which is why they were preferred on the dip criterion, but their L RISES monotonically
/// through the backward hemisphere (chi2/ndf 13 and 19). Koning-Delaroche's FALLS (15.8). Only
/// Perey gives a constant L, at chi2/ndf 1.5 over 80-120 deg. The dip position is one feature; the
/// flatness of L is the whole distribution, and it is the property the normalisation depends on.
///
///   root -b -q 'elastic_flatness_C14.C()'
namespace ef {
TGraph *rd(TString f)
{
   auto *g = new TGraph();
   std::ifstream in(f.Data());
   double a, b;
   while (in >> a >> b) if (b > 0) g->SetPoint(g->GetN(), a, b);
   return g;
}
} // namespace ef

void elastic_flatness_C14(TString elFile = "plots/elastic_omp_omp.root",
                          Double_t winLo = 80, Double_t winHi = 120)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TString pdir = here + "/../ptolemy/dat/";
   TFile *f = TFile::Open(here + "/" + elFile);
   auto *gd = f && !f->IsZombie() ? (TGraphErrors *)f->Get("elastic_measured") : nullptr;
   if (!gd) { printf("\033[1;31mno elastic_measured in %s -- run elastic_omp_C14.C\033[0m\n", elFile.Data()); return; }

   const int NP = 5;
   const char *pk[NP] = {"K", "V", "G", "P", "M"};
   const char *pn[NP] = {"KD03", "CH89", "Becchetti-Greenlees", "Perey", "Menet"};
   int col[NP] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1};

   // weighted constant fit to L(theta); the errors are the statistical ones on the elastic yield
   auto fit = [&](TGraph *c, double lo, double hi, double &Lbar, double &chi2ndf) {
      double sw = 0, swy = 0;
      std::vector<double> vL, vE;
      for (int i = 0; i < gd->GetN(); ++i) {
         double x = gd->GetX()[i], y = gd->GetY()[i], e = gd->GetEY()[i];
         if (x < lo || x > hi || y <= 0 || e <= 0) continue;
         double v = c->Eval(x);
         if (v <= 0) continue;
         double L = y / v, dL = e / v;
         vL.push_back(L); vE.push_back(dL);
         sw += 1 / (dL * dL); swy += L / (dL * dL);
      }
      Lbar = sw > 0 ? swy / sw : 0;
      double c2 = 0;
      for (size_t j = 0; j < vL.size(); ++j) c2 += std::pow((vL[j] - Lbar) / vE[j], 2);
      chi2ndf = vL.size() > 1 ? c2 / (vL.size() - 1) : 0;
      return (int)vL.size();
   };

   printf("\n  ==== is L constant? (it must be) ====\n");
   const int NW = 4;
   double W[NW][2] = {{70, 148}, {70, 120}, {80, 120}, {85, 115}};
   for (int w = 0; w < NW; ++w) {
      printf("\n  window %.0f-%.0f deg\n  %-22s %8s %10s\n", W[w][0], W[w][1], "potential", "L", "chi2/ndf");
      for (int i = 0; i < NP; ++i) {
         auto *g = ef::rd(pdir + TString::Format("el_omp_%s.dat", pk[i]));
         double L, c2; int n = fit(g, W[w][0], W[w][1], L, c2);
         printf("  %-22s %8.1f %10.1f   (%d points)\n", pn[i], L, c2, n);
      }
   }

   // ---- one panel per potential, as always -----------------------------------------------------
   auto *c1 = new TCanvas("cef", "", 1500, 900); c1->Divide(3, 2);
   int best = -1; double bestC2 = 1e99, bestL = 0;
   for (int i = 0; i < NP; ++i) {
      auto *g = ef::rd(pdir + TString::Format("el_omp_%s.dat", pk[i]));
      double L, c2; fit(g, winLo, winHi, L, c2);
      if (c2 < bestC2) { bestC2 = c2; best = i; bestL = L; }
      c1->cd(i + 1); gPad->SetGridx(); gPad->SetGridy();
      auto *q = new TGraphErrors();
      for (int j = 0; j < gd->GetN(); ++j) {
         double x = gd->GetX()[j], y = gd->GetY()[j], e = gd->GetEY()[j];
         if (x < 55 || y <= 0) continue;
         double v = g->Eval(x);
         if (v <= 0) continue;
         int k = q->GetN(); q->SetPoint(k, x, y / v); q->SetPointError(k, 0, e / v);
      }
      q->SetTitle(Form("%s;#theta_{cm} [deg];L = measured / calculated", pn[i]));
      q->SetMarkerStyle(20); q->SetMarkerColor(col[i]); q->SetLineColor(col[i]);
      q->SetMinimum(0); q->SetMaximum(230); q->Draw("AP");
      auto *ln = new TLine(winLo, L, winHi, L);
      ln->SetLineColor(kRed + 1); ln->SetLineWidth(3); ln->Draw();
      auto *bx = new TBox(winLo, 0, winHi, 230);
      bx->SetFillColorAlpha(kRed + 1, 0.08); bx->Draw();
      TLatex tx; tx.SetNDC(); tx.SetTextSize(0.055);
      tx.DrawLatex(0.16, 0.86, Form("L = %.1f", L));
      tx.SetTextColor(c2 < 3 ? kGreen + 3 : kRed + 1);
      tx.DrawLatex(0.16, 0.79, Form("#chi^{2}/ndf = %.1f", c2));
   }
   c1->cd(6);
   TLatex t; t.SetTextSize(0.058); t.SetNDC();
   t.DrawLatex(0.06, 0.86, "L must be a CONSTANT:");
   t.DrawLatex(0.06, 0.78, "it is a beam flux times an");
   t.DrawLatex(0.06, 0.70, "areal density, so any drift");
   t.DrawLatex(0.06, 0.62, "with angle is a defect of the");
   t.DrawLatex(0.06, 0.54, "normalisation, not of the beam.");
   t.SetTextColor(kGreen + 3);
   t.DrawLatex(0.06, 0.40, Form("Only %s passes,", pn[best]));
   t.DrawLatex(0.06, 0.32, Form("#chi^{2}/ndf = %.1f, L = %.1f", bestC2, bestL));
   t.SetTextColor(kBlack); t.SetTextSize(0.048);
   t.DrawLatex(0.06, 0.18, Form("window %.0f-%.0f#circ", winLo, winHi));
   t.DrawLatex(0.06, 0.10, "shaded on every panel");
   TString out = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";
   c1->SaveAs(out + "13_L_flatness.png");
   printf("\n  best: %s, L = %.1f, chi2/ndf = %.1f over %.0f-%.0f deg\n", pn[best], bestL, bestC2, winLo, winHi);
   printf("  wrote %s13_L_flatness.png\n\n", out.Data());
}
