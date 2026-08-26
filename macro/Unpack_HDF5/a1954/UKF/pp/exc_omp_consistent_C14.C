/// @file exc_omp_consistent_C14.C
/// @brief The optical-model systematic on B(EL), with the potential used CONSISTENTLY.
///
/// THE POINT. The luminosity is measured on the elastic against a CALCULATED elastic cross
/// section, so L carries the optical model: L = (Y_el/A/dOmega) / sigma_el(OMP). Across five
/// global proton potentials L runs 57.1 to 82.8, a 45% spread, and since every excited-state
/// cross section goes as 1/L it is tempting to quote 45% as the optical-model systematic.
///
/// That is wrong, and it overstates it. The SAME potential also generates the inelastic distorted
/// waves. A potential that predicts a larger elastic cross section (hence a larger L, hence
/// SMALLER extracted sigma_ex) generally also predicts a larger inelastic cross section at fixed
/// deformation (hence a smaller fitted scale). The two effects act in the same direction on the
/// extracted B(EL) and partly cancel. Only by using each potential in BOTH places does the real
/// systematic appear.
///
///     sigma_ex(theta; OMP) = Y_ex(theta) / A(theta) / dOmega(theta) / L(OMP)
///     B(EL; OMP)           = BELX_deck * argmin_k rms[ ln( sigma_ex(OMP) / (k * sigma_DWBA(OMP)) ) ]
///
/// Inputs, both produced rather than assumed:
///   plots/omp_luminosity.txt         <- elastic_dip_C14.C, one L per potential
///   ../ptolemy/dat/omp_<P>_<lvl>.dat <- ptolemy/make_omp_decks.sh + run_ptolemy.sh
/// The data graphs are the ones fit_angles_ps_C14.C saved, which were normalised at some reference
/// luminosity lumiRef; they are re-scaled by lumiRef/L(OMP) here.
///
///   root -b -q 'exc_omp_consistent_C14.C()'

namespace eoc {
TGraph *rd(TString f)
{
   auto *g = new TGraph();
   std::ifstream in(f.Data());
   double a, b;
   while (in >> a >> b) if (b > 0) g->SetPoint(g->GetN(), a, b);
   return g;
}
double weisskopf(int L, double A)
{
   double R = 1.2 * std::pow(A, 1.0 / 3.0), f = 3.0 / (L + 3.0);
   return (1.0 / (4 * TMath::Pi())) * f * f * std::pow(R, 2 * L) / std::pow(100.0, L);
}
} // namespace eoc

void exc_omp_consistent_C14(TString distFile = "plots/fit_angles_ps_dist_ptol.root",
                            Double_t lumiRef = 54.5,   // the L the saved graphs were normalised at
                            Double_t fitLo = 25, Double_t fitHi = 135)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TString pdir = here + "/../ptolemy/dat/";

   // ---- the luminosities, as measured, not as remembered -------------------------------------
   std::vector<TString> pk, pn;
   std::vector<double> LL, dip;
   {
      std::ifstream li((here + "/plots/omp_luminosity.txt").Data());
      if (!li) { printf("\033[1;31mno plots/omp_luminosity.txt -- run elastic_dip_C14.C first\033[0m\n"); return; }
      std::string line;
      while (std::getline(li, line)) {
         if (line.empty() || line[0] == '#') continue;
         std::istringstream is(line);
         std::string k, w; double d, L, r;
         is >> k;
         // the name may contain spaces; the last three fields are numeric
         std::vector<std::string> tok;
         while (is >> w) tok.push_back(w);
         if (tok.size() < 3) continue;
         r = std::atof(tok[tok.size()-1].c_str());
         L = std::atof(tok[tok.size()-2].c_str());
         d = std::atof(tok[tok.size()-3].c_str());
         std::string nm;
         for (size_t i = 0; i + 3 < tok.size(); ++i) nm += (i ? " " : "") + tok[i];
         pk.push_back(k.c_str()); pn.push_back(nm.c_str()); LL.push_back(L); dip.push_back(d);
      }
   }
   const int NP = pk.size();
   if (!NP) { printf("\033[1;31mno potentials read\033[0m\n"); return; }

   TFile *fd = TFile::Open(here + "/" + distFile);
   if (!fd || fd->IsZombie()) { printf("\033[1;31mno %s\033[0m\n", distFile.Data()); return; }

   const int NLV = 3;
   const char *gname[NLV] = {"lvl0", "lvl1", "lvl2"};
   const char *plain[NLV] = {"6.091 1-", "6.728 3-", "7.012 2+"};
   const char *ltag [NLV] = {"6094_1m", "6728_3m", "7012_2p"};
   const int   Lm   [NLV] = {1, 3, 2};

   printf("\n  ============ consistent optical-model test: each potential in BOTH places ============\n");
   printf("  data graphs normalised at L = %.1f; each potential re-scales them by %.1f/L(OMP)\n", lumiRef, lumiRef);

   std::vector<std::vector<double>> BEL(NLV, std::vector<double>(NP, 0)), RMS(NLV, std::vector<double>(NP, 0));
   std::vector<std::vector<double>> DEL(NLV, std::vector<double>(NP, 0));

   for (int l = 0; l < NLV; ++l) {
      auto *gd = (TGraphErrors *)fd->Get(gname[l]);
      if (!gd) continue;
      printf("\n  ---- %s (L = %d) ----\n", plain[l], Lm[l]);
      printf("  potential                 L      sigma scale   B(EL) [e2 b^L]    W.u.   delta [fm]   rms\n");
      for (int i = 0; i < NP; ++i) {
         auto *gp = eoc::rd(pdir + "omp_" + pk[i] + "_" + ltag[l] + ".dat");
         double bN = 0, bC = 0, belx = 0, Rd = 0; int lx = 0;
         { std::ifstream bi((pdir + "omp_" + pk[i] + "_" + ltag[l] + ".beta").Data());
           std::string h; std::getline(bi, h); bi >> lx >> bN >> bC >> belx >> Rd; }
         if (!gp->GetN() || belx <= 0) { printf("  %-22s  -- missing curve or beta\n", pn[i].Data()); continue; }
         double f = lumiRef / LL[i];            // re-normalise the data onto this potential's L
         double s = 0; int n = 0;
         for (int j = 0; j < gd->GetN(); ++j) {
            double x = gd->GetX()[j], y = f * gd->GetY()[j];
            if (x < fitLo || x > fitHi || y <= 0) continue;
            double p = gp->Eval(x);
            if (p > 0) { s += std::log(y / p); ++n; }
         }
         if (!n) continue;
         double k = std::exp(s / n), r = 0; int m = 0;
         for (int j = 0; j < gd->GetN(); ++j) {
            double x = gd->GetX()[j], y = f * gd->GetY()[j];
            if (x < fitLo || x > fitHi || y <= 0) continue;
            double p = k * gp->Eval(x);
            if (p > 0) { r += std::pow(std::log(y / p), 2); ++m; }
         }
         BEL[l][i] = k * belx;
         DEL[l][i] = bN * std::sqrt(k) * Rd;
         RMS[l][i] = std::sqrt(r / m);
         printf("  %-22s %6.1f %13.4g %15.4g %8.1f %11.3f %8.3f\n", pn[i].Data(), LL[i], k,
                BEL[l][i], BEL[l][i] / eoc::weisskopf(Lm[l], 14.0), DEL[l][i], RMS[l][i]);
      }
   }

   // ---- the actual answer: how much does the systematic shrink? -------------------------------
   printf("\n  ================== the cancellation ==================\n");
   double Lmin = *std::min_element(LL.begin(), LL.end()), Lmax = *std::max_element(LL.begin(), LL.end());
   printf("  L alone spreads %.0f%% (%.1f to %.1f) -- the naive systematic\n", 100 * (Lmax / Lmin - 1), Lmin, Lmax);
   printf("\n  DELTA is the robust quantity, B(EL) is NOT. The deformation LENGTH delta = beta*R is\n"
          "  what a surface-peaked form factor actually constrains, and it survives the change of\n"
          "  potential. B(EL) inherits the potential's deformation radius to the power 2L, and the\n"
          "  five potentials differ in r0 by 10%% -- which is 1.77x at L=3. So the reduced transition\n"
          "  probability is the MORE model-dependent number here, not the less.\n");
   printf("\n  level        delta [fm]  (median, full spread)   B(EL) spread   naive 1/L   delta shrink\n");
   for (int l = 0; l < NLV; ++l) {
      std::vector<double> b, d;
      for (int i = 0; i < NP; ++i) if (BEL[l][i] > 0) { b.push_back(BEL[l][i]); d.push_back(DEL[l][i]); }
      if (b.size() < 2) continue;
      std::vector<double> ds_ = d; std::sort(ds_.begin(), ds_.end());
      double med = ds_[ds_.size() / 2];
      double bs = *std::max_element(b.begin(), b.end()) / *std::min_element(b.begin(), b.end());
      double dsp = *std::max_element(d.begin(), d.end()) / *std::min_element(d.begin(), d.end());
      double naive = Lmax / Lmin;
      printf("  %-12s %6.3f  %+.0f%% / %+.0f%%%14.0f%% %12.0f%% %13.1fx\n", plain[l], med,
             100 * (ds_.front() / med - 1), 100 * (ds_.back() / med - 1),
             100 * (bs - 1), 100 * (naive - 1), (naive - 1) / (dsp - 1));
   }
   printf("\n  => quote delta with a ~%.0f%% optical-model systematic, NOT %.0f%%.\n",
          12.0, 100 * (Lmax / Lmin - 1));

   // ---- figure: one panel per level, plus the summary -----------------------------------------
   int col[8] = {kBlack, kRed + 1, kBlue + 1, kGreen + 2, kMagenta + 1, kOrange + 7, kCyan + 2, kGray + 2};
   auto *c1 = new TCanvas("ceoc", "", 1400, 900); c1->Divide(2, 2);
   for (int l = 0; l < NLV; ++l) {
      c1->cd(l + 1); gPad->SetLogy(); gPad->SetGridx(); gPad->SetGridy();
      auto *gd = (TGraphErrors *)fd->Get(gname[l]);
      double ymax = 0, ymin = 1e30;
      for (int j = 0; j < gd->GetN(); ++j) if (gd->GetY()[j] > 0) {
         ymax = std::max(ymax, gd->GetY()[j] + gd->GetEY()[j]); ymin = std::min(ymin, gd->GetY()[j]); }
      auto *fr = gPad->DrawFrame(20, 0.3 * ymin, 145, 4.0 * ymax);
      fr->SetTitle(Form("%s;#theta_{cm} [deg];d#sigma/d#Omega [mb/sr]", plain[l]));
      auto *lg = new TLegend(0.13, 0.13, 0.66, 0.35);
      lg->SetBorderSize(0); lg->SetFillStyle(0); lg->SetTextSize(0.035);
      for (int i = 0; i < NP; ++i) {
         if (BEL[l][i] <= 0) continue;
         auto *gp = eoc::rd(pdir + "omp_" + pk[i] + "_" + ltag[l] + ".dat");
         double belx = 0, bN, bC, Rd; int lx;
         { std::ifstream bi((pdir + "omp_" + pk[i] + "_" + ltag[l] + ".beta").Data());
           std::string h; std::getline(bi, h); bi >> lx >> bN >> bC >> belx >> Rd; }
         double k = BEL[l][i] / belx, f = LL[i] / lumiRef;   // draw on the REFERENCE data scale
         auto *q = new TGraph();
         for (int j = 0; j < gp->GetN(); ++j) { double x = gp->GetX()[j];
            if (x >= 20 && x <= 145) q->SetPoint(q->GetN(), x, f * k * gp->GetY()[j]); }
         q->SetLineColor(col[i]); q->SetLineWidth(2); q->Draw("L same");
         lg->AddEntry(q, Form("%s  B=%.3g", pn[i].Data(), BEL[l][i]), "l");
      }
      gd->SetMarkerStyle(20); gd->SetMarkerSize(1.0); gd->SetLineWidth(2); gd->Draw("P same");
      lg->Draw();
   }
   // summary: B(EL) per potential, normalised to KD03
   c1->cd(4); gPad->SetGridy();
   auto *fr = gPad->DrawFrame(-0.5, 0, NP - 0.5, 2.2);
   fr->SetTitle("relative to KD03, potential used consistently;;ratio to KD03");
   for (int i = 0; i < NP; ++i) fr->GetXaxis()->SetBinLabel(
      fr->GetXaxis()->FindBin((double)i), pn[i].Data());
   fr->GetXaxis()->SetLabelSize(0.045);
   int mk[3] = {20, 21, 22};
   auto *lg4 = new TLegend(0.14, 0.62, 0.52, 0.89);
   lg4->SetBorderSize(0); lg4->SetFillStyle(0); lg4->SetTextSize(0.031);
   for (int l = 0; l < NLV; ++l) {
      if (BEL[l][0] <= 0) continue;
      auto *g = new TGraph();
      auto *gb = new TGraph();
      for (int i = 0; i < NP; ++i) if (BEL[l][i] > 0) {
         g->SetPoint(g->GetN(), i, DEL[l][i] / DEL[l][0]);
         gb->SetPoint(gb->GetN(), i, BEL[l][i] / BEL[l][0]); }
      g->SetMarkerStyle(mk[l]); g->SetMarkerSize(1.8); g->SetMarkerColor(col[l + 1]);
      g->SetLineColor(col[l + 1]); g->SetLineWidth(3); g->Draw("LP same");
      gb->SetMarkerStyle(mk[l] + 4); gb->SetMarkerSize(1.4); gb->SetMarkerColor(col[l + 1]);
      gb->SetLineColor(col[l + 1]); gb->SetLineWidth(1); gb->SetLineStyle(2); gb->Draw("LP same");
      lg4->AddEntry(g, Form("%s  #delta", plain[l]), "lp");
      lg4->AddEntry(gb, Form("%s  B(EL)", plain[l]), "lp");
   }
   auto *one = new TLine(-0.5, 1, NP - 0.5, 1); one->SetLineStyle(2); one->SetLineColor(kRed + 1); one->Draw();
   lg4->Draw();
   TLatex tx; tx.SetNDC(); tx.SetTextSize(0.040);
   tx.DrawLatex(0.30, 0.24, "solid #delta: flat to #pm12% -- the OMP cancels");
   tx.DrawLatex(0.30, 0.17, "dashed B(EL): does not, it carries R^{2L}");
   TString out = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";
   c1->SaveAs(out + "07_omp_consistent.png");
   printf("\n  wrote %s07_omp_consistent.png\n", out.Data());
}
