/// @file ch13_figs_C14.C
/// @brief Regenerate the Chapter 13 figures and table on the adopted normalisation (Perey).
///
/// WHY NOT Ldep_C14.C AND absnorm_C14.C. Those read the 2026-08 genfit production
/// (acc_gf_z10_400, proton_kin_300gfx_nc) and a hardcoded Koning--Delaroche elastic. Re-running
/// them with a different potential would mix that production with the CATIMA one the rest of the
/// analysis now uses -- two different acceptances over two different vertex slabs -- which is the
/// same class of silent inconsistency that produced the original factor 5-8. These figures are
/// therefore built from the SAME chain that produced the numbers in Section 13.5 and in the
/// companion document, so the three cannot disagree.
///
/// Two figures:
///   ldep_perey    L(theta) = measured / calculated for the adopted potential, with the
///                 80-120 deg window. This replaces Ldep_C14.png, whose 40-60 deg band sits on
///                 the shoulder of a diffraction minimum Koning--Delaroche misplaces by 5 deg.
///   absnorm_perey the three multiplet levels in absolute mb/sr against their own Perey DWBA
///                 curves, each on one free scale which IS the deformation.
///
/// and the numbers for the replacement of tab:absnorm.
///
///   root -b -q 'ch13_figs_C14.C()'
namespace c13 {
TGraph *rd(TString f)
{
   auto *g = new TGraph();
   std::ifstream in(f.Data());
   double a, b;
   while (in >> a >> b) if (b > 0) g->SetPoint(g->GetN(), a, b);
   return g;
}
} // namespace c13

void ch13_figs_C14(TString elFile = "plots/elastic_omp_omp.root",
                   TString distFile = "plots/fit_angles_ps_dist_pr00.root",
                   TString omp = "P", TString ompName = "Perey",
                   Double_t winLo = 80, Double_t winHi = 120,
                   Double_t fitLo = 25, Double_t fitHi = 135)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TString pdir = here + "/../ptolemy/dat/";
   TString out = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";

   // ---------- figure 1: L(theta) ----------
   TFile *fe = TFile::Open(here + "/" + elFile);
   auto *gd = fe && !fe->IsZombie() ? (TGraphErrors *)fe->Get("elastic_measured") : nullptr;
   auto *gc = c13::rd(pdir + "el_omp_" + omp + ".dat");
   if (!gd || !gc->GetN()) { printf("\033[1;31mmissing elastic or curve\033[0m\n"); return; }

   double sw = 0, swy = 0;
   std::vector<double> vL, vE, vT;
   for (int i = 0; i < gd->GetN(); ++i) {
      double x = gd->GetX()[i], y = gd->GetY()[i], e = gd->GetEY()[i];
      if (y <= 0 || e <= 0) continue;
      double v = gc->Eval(x);
      if (v <= 0) continue;
      vT.push_back(x); vL.push_back(y / v); vE.push_back(e / v);
      if (x >= winLo && x <= winHi) { sw += 1 / std::pow(e / v, 2); swy += (y / v) / std::pow(e / v, 2); }
   }
   double Lbar = sw > 0 ? swy / sw : 0, c2 = 0; int nw = 0;
   for (size_t j = 0; j < vT.size(); ++j)
      if (vT[j] >= winLo && vT[j] <= winHi) { c2 += std::pow((vL[j] - Lbar) / vE[j], 2); ++nw; }
   double c2n = nw > 1 ? c2 / (nw - 1) : 0;

   auto *c1 = new TCanvas("c13a", "", 1300, 560); c1->Divide(2, 1);
   c1->cd(1); gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
   auto *fr = gPad->DrawFrame(15, 5, 152, 1e5);
   fr->SetTitle(Form("elastic, measured against %s;#theta_{cm} [deg];d#sigma/d#Omega [arb.]", ompName.Data()));
   auto *q = new TGraph();
   for (int j = 0; j < gc->GetN(); ++j) { double x = gc->GetX()[j];
      if (x >= 15 && x <= 152) q->SetPoint(q->GetN(), x, Lbar * gc->GetY()[j]); }
   q->SetLineColor(kRed + 1); q->SetLineWidth(3); q->Draw("L same");
   gd->SetMarkerStyle(20); gd->SetMarkerSize(1.0); gd->Draw("P same");
   c1->cd(2); gPad->SetGridx(); gPad->SetGridy();
   auto *gl = new TGraphErrors(vT.size(), &vT[0], &vL[0], nullptr, &vE[0]);
   gl->SetTitle("L = measured / calculated, which must be a constant;#theta_{cm} [deg];L [counts/mb]");
   gl->SetMarkerStyle(20); gl->SetMinimum(0); gl->SetMaximum(250); gl->Draw("AP");
   auto *bx = new TBox(winLo, 0, winHi, 250); bx->SetFillColorAlpha(kRed + 1, 0.10); bx->Draw();
   auto *ln = new TLine(winLo, Lbar, winHi, Lbar); ln->SetLineColor(kRed + 1); ln->SetLineWidth(3); ln->Draw();
   gl->Draw("P same");
   TLatex tx; tx.SetNDC(); tx.SetTextSize(0.042);
   tx.DrawLatex(0.14, 0.86, Form("L = %.1f counts/mb over %.0f-%.0f#circ", Lbar, winLo, winHi));
   tx.DrawLatex(0.14, 0.80, Form("#chi^{2}/ndf = %.1f for a constant", c2n));
   tx.SetTextColor(kGray + 2); tx.SetTextSize(0.036);
   tx.DrawLatex(0.14, 0.72, "the old 40-60#circ band sits on the");
   tx.DrawLatex(0.14, 0.67, "shoulder of the diffraction minimum");
   c1->SaveAs(out + "14_ldep_" + omp + ".png");

   // ---------- figure 2 + the table ----------
   TFile *fd = TFile::Open(here + "/" + distFile);
   if (!fd || fd->IsZombie()) { printf("\033[1;31mno %s\033[0m\n", distFile.Data()); return; }
   const int NL = 3;
   const char *gn[NL] = {"lvl0", "lvl1", "lvl2"};
   const char *nm[NL] = {"6.091 MeV 1^{-}", "6.728 MeV 3^{-}", "7.012 MeV 2^{+}"};
   const char *pl[NL] = {"6.091 1-", "6.728 3-", "7.012 2+"};
   const char *pf[NL] = {"6094_1m", "6728_3m", "7012_2p"};
   const int   Lm[NL] = {1, 3, 2};
   int col[NL] = {kAzure + 2, kGreen + 3, kOrange + 7};

   printf("\n  ==== replacement for tab:absnorm (%s, L = %.1f) ====\n", ompName.Data(), Lbar);
   printf("  %-10s %3s %14s %10s %8s %8s\n", "level", "L", "sigma(25-135)", "delta[fm]", "beta", "rms");
   auto *c2c = new TCanvas("c13b", "", 1400, 480); c2c->Divide(3, 1);
   for (int i = 0; i < NL; ++i) {
      auto *g = (TGraphErrors *)fd->Get(gn[i]);
      auto *gp = c13::rd(pdir + "omp_" + omp + "_" + pf[i] + ".dat");
      double bN, bC, belx, R; int lx;
      { std::ifstream b((pdir + "omp_" + omp + "_" + pf[i] + ".beta").Data());
        std::string h; std::getline(b, h); b >> lx >> bN >> bC >> belx >> R; }
      if (!g || !gp->GetN() || belx <= 0) { printf("  %-10s missing\n", pl[i]); continue; }
      double s = 0, sig = 0; int n = 0;
      for (int j = 0; j < g->GetN(); ++j) {
         double x = g->GetX()[j], y = g->GetY()[j];
         if (x < fitLo || x > fitHi || y <= 0) continue;
         double p = gp->Eval(x);
         if (p > 0) { s += std::log(y / p); ++n; }
         sig += y * 2 * TMath::Pi() * std::sin(x * TMath::DegToRad()) * 10 * TMath::DegToRad();
      }
      double k = std::exp(s / n), r = 0; int m = 0;
      for (int j = 0; j < g->GetN(); ++j) {
         double x = g->GetX()[j], y = g->GetY()[j];
         if (x < fitLo || x > fitHi || y <= 0) continue;
         double p = k * gp->Eval(x);
         if (p > 0) { r += std::pow(std::log(y / p), 2); ++m; }
      }
      double beta = bN * std::sqrt(k), delta = beta * R;
      printf("  %-10s %3d %14.1f %10.3f %8.3f %8.3f\n", pl[i], Lm[i], sig, delta, beta, std::sqrt(r / m));

      c2c->cd(i + 1); gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
      // DROP THE RAILED POINTS. Where a level is too weak to separate -- the 6.091 and the 6.728
      // beyond ~105 deg in the 5 deg binning -- the amplitude is pushed onto its lower limit and
      // comes back at ~1e-10 rather than exactly zero, so it survives the y > 0 filter upstream.
      // Those points carry no information, they are already excluded from the normalisation by the
      // fit range, and on a log axis they drag the frame down ten decades and flatten every real
      // point into the top sliver. A point is called railed if it is below 1e-3 of the largest
      // measured value for that level; the count is reported in the panel title so the removal is
      // visible rather than silent.
      double ypk = 0;
      for (int j = 0; j < g->GetN(); ++j) ypk = std::max(ypk, g->GetY()[j]);
      int nrail = 0;
      for (int j = g->GetN() - 1; j >= 0; --j)
         if (g->GetY()[j] <= 0 || g->GetY()[j] < 1e-3 * ypk) { g->RemovePoint(j); ++nrail; }
      double ymax = 0, ymin = 1e30;
      for (int j = 0; j < g->GetN(); ++j) if (g->GetY()[j] > 0) {
         ymax = std::max(ymax, g->GetY()[j] + g->GetEY()[j]); ymin = std::min(ymin, g->GetY()[j]); }
      if (ymin > ymax) { ymin = 0.1; ymax = 1.0; }
      auto *f2 = gPad->DrawFrame(20, 0.4 * ymin, 145, 4.0 * ymax);
      f2->SetTitle(nrail ? Form("%s   [%d railed bin%s not shown];#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]",
                                nm[i], nrail, nrail > 1 ? "s" : "")
                         : Form("%s;#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]", nm[i]));
      auto *qq = new TGraph();
      for (int j = 0; j < gp->GetN(); ++j) { double x = gp->GetX()[j];
         if (x >= 20 && x <= 145) qq->SetPoint(qq->GetN(), x, k * gp->GetY()[j]); }
      qq->SetLineColor(col[i]); qq->SetLineWidth(3); qq->Draw("L same");
      g->SetMarkerStyle(20); g->SetMarkerSize(1.1); g->SetLineWidth(2); g->Draw("P same");
      TLatex t2; t2.SetNDC(); t2.SetTextSize(0.048);
      t2.DrawLatex(0.16, 0.26, Form("#delta = %.3f fm", delta));
      t2.DrawLatex(0.16, 0.19, Form("#beta = %.3f", beta));
   }
   c2c->SaveAs(out + "15_absnorm_" + omp + ".png");
   printf("\n  wrote %s14_ldep_%s.png and %s15_absnorm_%s.png\n\n", out.Data(), omp.Data(), out.Data(), omp.Data());
}
