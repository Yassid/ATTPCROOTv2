/// @file exc_ptolemy_C14.C
/// @brief The 14C excited states against PtolemyCpp, and B(EL) instead of a deformation length.
///
/// What this adds over exc_vs_fresco_C14.C. FRESCO's collective coupling is parametrised by a
/// deformation LENGTH, so a fitted scale k only gives delta = delta_ref * sqrt(k) -- a number that
/// depends on the reference chosen and is awkward to compare with anything. Ptolemy's INELOCA1
/// coupling is parametrised by BELX = B(EL) itself, in e^2 barn^L, and the cross section is
/// LINEAR in it (verified to 3e-5 by rerunning at BELX=0.01, see ptolemy/inputs/*_linearity.in).
/// So the same fitted scale gives the reduced transition probability directly:
///
///     B(EL) = k * BELX_input
///
/// which can be quoted in Weisskopf units and compared with the literature. Ptolemy also prints
/// the nuclear deformation and radius it derived, so delta = beta_N * R comes out too and can be
/// checked against the FRESCO number.
///
/// CAVEATS, both real:
///  - The collective model is for even-even spherical HEAVY nuclei. 14C at 6-7 MeV is neither, so
///    these B(EL) are "what a collective model would need", not measured collectivity.
///  - 7.341 is 2-, unnatural parity. A one-step collective calculation CANNOT populate it: Ptolemy
///    returns no cross section at all for a 2- state, which is the honest answer. L=1 and L=3 are
///    run as SHAPE proxies only and their B(EL) is meaningless.
///
///   root -b -q 'exc_ptolemy_C14.C()'

namespace ep {
TGraph *rd(TString f)
{
   auto *g = new TGraph();
   std::ifstream in(f.Data());
   double a, b;
   while (in >> a >> b) if (b > 0) g->SetPoint(g->GetN(), a, b);
   if (!g->GetN()) printf("\033[1;31mcannot read %s\033[0m\n", f.Data());
   return g;
}
// Weisskopf single-particle estimate, B(EL)_W = (1/4pi)[3/(L+3)]^2 e^2 R^(2L), R = 1.2 A^(1/3) fm.
// Returned in e^2 barn^L to match Ptolemy's BELX, i.e. divided by 100^L.
double weisskopf(int L, double A)
{
   double R = 1.2 * std::pow(A, 1.0 / 3.0), f = 3.0 / (L + 3.0);
   return (1.0 / (4 * TMath::Pi())) * f * f * std::pow(R, 2 * L) / std::pow(100.0, L);
}
} // namespace ep

void exc_ptolemy_C14(TString distFile = "plots/fit_angles_ps_dist_ptol.root",
                     Double_t fitLo = 25, Double_t fitHi = 135)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TString pdir = here + "/../ptolemy/dat/", fdir = here + "/../fresco/outputs/";
   TFile *fd = TFile::Open(here + "/" + distFile);
   if (!fd || fd->IsZombie()) { printf("\033[1;31mno %s -- run fit_angles_ps_C14 first\033[0m\n", distFile.Data()); return; }
   if (auto *pv = (TNamed *)fd->Get("provenance")) printf("\n  data: %s\n", pv->GetTitle());

   const int NL = 4;
   const char *gname[NL] = {"lvl0", "lvl1", "lvl2", "lvl3"};
   const char *label[NL] = {"6.091 1^{-}", "6.728 3^{-}", "7.012 2^{+}", "7.341 2^{-}"};
   const char *plain[NL] = {"6.091 1-", "6.728 3-", "7.012 2+", "7.341 2-"};
   const char *pfile[NL] = {"exc_6094_1m", "exc_6728_3m", "exc_7012_2p", "exc_7341_1m"};
   const char *ffile[NL] = {"p14C_inel_161_6094_L1_dsdo_ex2.dat", "p14C_inel_161_6728_L3_dsdo_ex2.dat",
                            "p14C_inel_161_7012_L2_dsdo_ex2.dat", "p14C_inel_161_7341_L1_dsdo_ex2.dat"};
   const int Lmult[NL] = {1, 3, 2, 1};
   const bool proxy[NL] = {false, false, false, true};   // 7.341 is unnatural parity
   // Ptolemy's nuclear deformation, Coulomb deformation and real deformation radius at the BELX
   // each deck was run at, carried into dat/*.beta by parse_ptolemy.py. delta = beta_N * R, and
   // both scale as sqrt(k). These MUST be read, not assumed -- beta_N differs by a factor 8
   // between the L=1 and L=3 couplings at the same BELX.
   double betaN[NL] = {0}, betaC[NL] = {0}, Rdef[NL] = {0}, belx[NL] = {0};
   for (int i = 0; i < NL; ++i) {
      std::ifstream bi((pdir + TString(pfile[i]) + ".beta").Data());
      std::string hdr; std::getline(bi, hdr);
      int lx = 0;
      if (!(bi >> lx >> betaN[i] >> betaC[i] >> belx[i] >> Rdef[i])) {
         printf("\033[1;31m  no %s.beta -- rerun ptolemy/run_ptolemy.sh\033[0m\n", pfile[i]);
         return;
      }
      if (lx != Lmult[i]) printf("\033[1;33m  %s: deck multipole LX=%d, expected %d\033[0m\n",
                                 pfile[i], lx, Lmult[i]);
   }

   printf("\n  ================= 14C(p,p') excited states vs PtolemyCpp (INELOCA1, KD03) =================\n");
   printf("  level         L   B(EL) [e2 b^L]   in W.u.    beta_N   delta [fm]   rms ln(d/f)   FRESCO delta\n");

   auto *c1 = new TCanvas("cep", "", 1400, 900); c1->Divide(2, 2);
   std::vector<double> RMSP(NL, 0), RMSF(NL, 0);
   for (int i = 0; i < NL; ++i) {
      auto *gd = (TGraphErrors *)fd->Get(gname[i]);
      auto *gp = ep::rd(pdir + TString(pfile[i]) + ".dat");
      auto *gf = ep::rd(fdir + ffile[i]);
      if (!gd || !gp->GetN()) continue;

      // log-space scale on both codes, over the same angles
      double sp = 0, sf = 0; int n = 0;
      for (int j = 0; j < gd->GetN(); ++j) {
         double x = gd->GetX()[j], y = gd->GetY()[j];
         if (x < fitLo || x > fitHi || y <= 0) continue;
         double p = gp->Eval(x), f = gf->GetN() ? gf->Eval(x) : 0;
         if (p <= 0) continue;
         sp += std::log(y / p);
         if (f > 0) sf += std::log(y / f);
         ++n;
      }
      double kp = std::exp(sp / n), kf = std::exp(sf / n);
      double rp = 0, rf = 0; int m = 0;
      for (int j = 0; j < gd->GetN(); ++j) {
         double x = gd->GetX()[j], y = gd->GetY()[j];
         if (x < fitLo || x > fitHi || y <= 0) continue;
         double p = kp * gp->Eval(x), f = kf * (gf->GetN() ? gf->Eval(x) : 0);
         if (p <= 0) continue;
         rp += std::pow(std::log(y / p), 2);
         if (f > 0) rf += std::pow(std::log(y / f), 2);
         ++m;
      }
      RMSP[i] = std::sqrt(rp / m); RMSF[i] = std::sqrt(rf / m);

      double BEL = kp * belx[i], wu = BEL / ep::weisskopf(Lmult[i], 14.0);
      double bN = betaN[i] * std::sqrt(kp), del = bN * Rdef[i];
      double dFR = 0.281 * std::sqrt(kf);
      if (proxy[i])
         printf("  %-12s (%d)  --- shape proxy only, 2- cannot be reached in one step ---  %11.3f %12.3f\n",
                plain[i], Lmult[i], RMSP[i], dFR);
      else
         printf("  %-12s  %d  %14.4g %10.1f %8.3f %11.3f %13.3f %12.3f\n",
                plain[i], Lmult[i], BEL, wu, bN, del, RMSP[i], dFR);

      c1->cd(i + 1); gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
      double ymax = 0, ymin = 1e30;
      for (int j = 0; j < gd->GetN(); ++j) if (gd->GetY()[j] > 0) {
         ymax = std::max(ymax, gd->GetY()[j] + gd->GetEY()[j]); ymin = std::min(ymin, gd->GetY()[j]); }
      auto *fr = gPad->DrawFrame(20, 0.4 * ymin, 145, 3.0 * ymax);
      fr->SetTitle(Form("%s%s;#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]",
                        label[i], proxy[i] ? "  (L=1 shape proxy)" : ""));
      auto *qp = new TGraph();
      for (int j = 0; j < gp->GetN(); ++j) { double x = gp->GetX()[j];
         if (x >= 20 && x <= 145) qp->SetPoint(qp->GetN(), x, kp * gp->GetY()[j]); }
      qp->SetLineColor(kRed + 1); qp->SetLineWidth(3); qp->Draw("L same");
      auto *qf = new TGraph();
      for (int j = 0; j < gf->GetN(); ++j) { double x = gf->GetX()[j];
         if (x >= 20 && x <= 145) qf->SetPoint(qf->GetN(), x, kf * gf->GetY()[j]); }
      qf->SetLineColor(kBlue + 1); qf->SetLineWidth(2); qf->SetLineStyle(2); qf->Draw("L same");
      gd->SetMarkerStyle(20); gd->SetMarkerSize(1.1); gd->SetLineWidth(2); gd->Draw("P same");
      auto *lg = new TLegend(0.14, 0.14, 0.62, 0.34);
      lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.040);
      lg->AddEntry(gd, "data", "pe");
      lg->AddEntry(qp, Form("Ptolemy INELOCA1  (rms %.3f)", RMSP[i]), "l");
      if (qf->GetN()) lg->AddEntry(qf, Form("FRESCO  (rms %.3f)", RMSF[i]), "l");
      lg->Draw();
      if (!proxy[i]) {
         TLatex tx; tx.SetNDC(); tx.SetTextSize(0.042);
         tx.DrawLatex(0.46, 0.86, Form("B(E%d) = %.3g e^{2}b^{%d}", Lmult[i], BEL, Lmult[i]));
         tx.DrawLatex(0.46, 0.80, Form("= %.0f W.u.", wu));
      }
   }
   TString out = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";
   c1->SaveAs(out + "06_excited_states_ptolemy.png");
   printf("\n  wrote %s06_excited_states_ptolemy.png\n", out.Data());
}
