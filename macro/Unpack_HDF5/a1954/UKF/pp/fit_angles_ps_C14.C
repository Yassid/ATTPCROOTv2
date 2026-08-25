/// @file fit_angles_ps_C14.C
/// @brief Per-angle version of the single 5-10.4 MeV fit: all 14C levels, the 14N (p,n) blend and
///        the 1n continuum in ONE model, per theta_cm bin, with every component drawn.
///
/// Range starts at 5.0 MeV ON PURPOSE. There are no 14C levels between the g.s. and 6.094, so
/// 5.0-6.0 is a level-free region that MEASURES the background. Fitting from 6.0 instead leaves the
/// background a free line extrapolated from where it cannot be measured, and it comes out steep
/// (13 counts/50 keV at 6 MeV against a true floor of 3.2) -- which then vacates the region above
/// Sn and lets the 1n continuum appear to be needed. From 5.0 the background is flat and low, and
/// the continuum goes to zero. Both are kept in the model here; which of them takes the yield above
/// Sn is a real degeneracy with the 14N blend and is NOT resolved by chi2.
///
/// STAGE 1 (angle integrated) fixes the shapes: positions at literature + the measured shift,
/// widths from sigma(Ex) = 0.132 + 0.0123 (Ex - 6.094), the 14N centroid/width at their measured
/// values. STAGE 2 floats ONLY the amplitudes per bin -- with a few hundred counts a bin cannot
/// re-measure a position or a width, and letting it try makes neighbouring peaks trade yield.
///
///   root -b -q 'fit_angles_ps_C14.C("plots/proton_kin_cat5_tc.root")'
namespace fa {
// HOW MANY LEVELS ARE ACTUALLY POPULATED IS A REAL QUESTION, NOT A DETAIL.
// The published analysis of these data (Ayyad et al.) uses FOUR components -- 6.09 (1-),
// 6.70 (3-), 7.00 (2+), 7.27 (2-). fit_states_C14.C added 6.589 (0+) and 6.903 (0-) because
// without them the fit has no way to describe the yield between the others. But per angle bin
// those two are NOT separable: at 100-120 deg the 6.903 took 169 +- 76 while 6.728 dropped to
// 73 +- 54, i.e. they swapped, with errors as large as the values.
// nLevels = 5 -> the published four plus 8.317   (recommended per angle)
// nLevels = 8 -> all seven plus 8.317            (angle integrated only)
int NL = 5;
double E[8]  = {6.091, 6.728, 7.012, 7.341, 8.317, 6.589, 6.903, 0};
double SG[8];
double shift = -0.011, muN = 9.178, sgN = 0.296;
TH1D *hps = nullptr;
double model(double *x, double *p)
{
   double s = p[NL + 2] + p[NL + 3] * x[0];
   for (int i = 0; i < NL; ++i) s += p[i] * std::exp(-0.5 * std::pow((x[0] - (E[i] + fa::shift)) / SG[i], 2));
   s += p[NL] * std::exp(-0.5 * std::pow((x[0] - muN) / sgN, 2));
   if (hps) s += p[NL + 1] * hps->Interpolate(x[0]);
   return s;
}
} // namespace fa

void fit_angles_ps_C14(TString cache = "plots/proton_kin_cat5_tc.root", Double_t chi2Cut = 5.0,
                       Double_t exLo = 5.0, Double_t exHi = 10.4, Double_t shift_ = -0.011,
                       Double_t sig0 = 0.132, Double_t dSig = 0.0123, Double_t muN_ = 9.178,
                       Double_t sgN_ = 0.296, Double_t cmMin = 20, Double_t cmMax = 140,
                       Double_t dcm = 20, Int_t minN = 40, TString tag = "cat5",
                       Int_t nLevels = 5, TString accDir = "/mnt/f/a1954_C14_acc_gf_nochi2/",
                       Bool_t doAcc = kTRUE)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TFile *fd = TFile::Open(here + "/" + cache);
   TNtuple *t = (TNtuple *)fd->Get("pk");
   TFile *fp = TFile::Open(here + "/plots/phasespace_1n_C14.root");
   TH1D *hps = fp && !fp->IsZombie() ? (TH1D *)fp->Get("hPS") : nullptr;
   if (hps) { hps->SetDirectory(nullptr); hps->Scale(1.0 / hps->Integral()); }

   using namespace fa;
   fa::NL = (nLevels >= 8) ? 7 : 4;      // 4 published + 8.317, or all 7 + 8.317
   // 8.317 always sits LAST among the fitted levels so its index is known for the acceptance
   if (fa::NL == 7) { double t7[7] = {6.091, 6.589, 6.728, 6.903, 7.012, 7.341, 8.317};
                      for (int i = 0; i < 7; ++i) fa::E[i] = t7[i]; }
   else             { double t5[5] = {6.091, 6.728, 7.012, 7.341, 8.317};
                      for (int i = 0; i < 5; ++i) fa::E[i] = t5[i]; fa::NL = 5; }
   fa::shift = shift_; fa::muN = muN_; fa::sgN = sgN_; fa::hps = hps;
   const int NL = fa::NL;
   std::vector<TString> jp;
   for (int i = 0; i < NL; ++i) jp.push_back(TString::Format("%.3f", fa::E[i]));
   const char *JPfull[8] = {"6.091 1-", "6.589 0+", "6.728 3-", "6.903 0-", "7.012 2+", "7.341 2-", "8.317 2+", ""};
   std::vector<TString> JPv;
   for (int i = 0; i < NL; ++i) {
      TString nm = TString::Format("%.3f", fa::E[i]);
      for (int k = 0; k < 7; ++k) if (TString(JPfull[k]).BeginsWith(nm)) nm = JPfull[k];
      JPv.push_back(nm);
   }
   std::vector<const char *> JP;
   for (auto &x : JPv) JP.push_back(x.Data());
   static const int COL[10] = {kBlue + 1, kAzure + 7, kGreen + 2, kOrange + 7, kMagenta + 1,
                               kViolet + 1, kRed + 1, kCyan + 2, kGray + 2, kBlack};
   for (int i = 0; i < NL; ++i) fa::SG[i] = sig0 + dSig * (fa::E[i] - 6.094);
   const int NPAR = NL + 4;   // levels, 14N, PS, bg0, bg1
   int nb = (int)std::lround((exHi - exLo) / 0.05);

   const int NB = (int)std::lround((cmMax - cmMin) / dcm);
   std::vector<std::vector<double>> Y(NL + 2), EY(NL + 2);   // + 14N, + PS
   std::vector<double> TH_;
   // Layout scaled to the number of bins. At 5 deg there are 24 of them and a 12x2 grid makes
   // every panel an unreadable sliver -- pick a near-square grid and size the canvas to it.
   int ncol = (int)std::ceil(std::sqrt((double)NB));
   int nrow = (int)std::ceil((double)NB / ncol);
   auto *cf = new TCanvas("cfits", "per-angle fits", 420 * ncol, 320 * nrow);
   cf->Divide(ncol, nrow, 0.002, 0.002);
   printf("\n  theta_cm |   N  |");
   for (int i = 0; i < NL; ++i) printf(" %10s |", JP[i]);
   printf("  14N  |  1n PS | chi2/ndf\n");

   for (int b = 0; b < NB; ++b) {
      double lo = cmMin + b * dcm, hi = lo + dcm;
      auto *h = new TH1D(Form("hb%d", b), Form("%.0f-%.0f deg;E_{x} [MeV];counts", lo, hi), nb, exLo, exHi);
      h->Sumw2();
      t->Draw(Form("ex>>hb%d", b), Form("chi2ndf<%g && thcm>=%g && thcm<%g", chi2Cut, lo, hi), "goff");
      if (h->Integral() < minN) { printf("  %3.0f-%3.0f  | %4.0f | too few\n", lo, hi, h->Integral()); continue; }
      auto *F = new TF1(Form("F%d", b), fa::model, exLo, exHi, NPAR);
      for (int i = 0; i < NL; ++i) {
         F->SetParameter(i, std::max(1.0, h->GetBinContent(h->FindBin(E[i] + fa::shift)) * 0.5));
         F->SetParLimits(i, 0, 1e5);
      }
      F->SetParameter(NL, 20);        F->SetParLimits(NL, 0, 1e5);
      F->SetParameter(NL + 1, h->Integral() * 0.05); F->SetParLimits(NL + 1, 0, 1e7);
      F->SetParameter(NL + 2, 2);     F->SetParameter(NL + 3, 0);
      TFitResultPtr r = h->Fit(F, "RQNSL");
      // RETRY BEFORE GIVING UP. At 25-30 deg the first pass returned chi2/ndf 0.78 with sane
      // VALUES but every error at ~12000 (and 590000 on the continuum) -- MIGRAD stopped somewhere
      // HESSE could not be computed. Refitting FROM that minimum converges: the errors come back
      // at 0.4-0.9 and the yields are unchanged. Dropping the bin was throwing away good data.
      for (int pass = 0; pass < 3; ++pass) {
         bool errBad = false;
         for (int i = 0; i < NL; ++i)
            if (F->GetParError(i) * SG[i] * std::sqrt(2 * TMath::Pi()) / h->GetBinWidth(1) > h->Integral())
               errBad = true;
         if (!errBad && r.Get() && r->CovMatrixStatus() >= 2) break;
         r = h->Fit(F, "RQNSL");          // restarts from the current parameters
      }
      // A CONVERGED FIT WITH A BROKEN COVARIANCE IS THE DANGEROUS CASE: at 25-30 deg every
      // component came back with an error of ~80000 on a yield of ~50, while chi2/ndf was 0.78 and
      // the VALUES were sane. One such bin sets the axis of every angular-distribution panel and
      // flattens all the real points to invisibility. Detected two ways: Minuit's own covariance
      // status, and an error larger than the total counts in the bin, which is impossible.
      int covStat = r.Get() ? r->CovMatrixStatus() : 0;
      bool badErr = false;
      for (int i = 0; i < NL; ++i)
         if (F->GetParError(i) * SG[i] * std::sqrt(2 * TMath::Pi()) / h->GetBinWidth(1) > h->Integral())
            badErr = true;
      bool bad = (covStat < 2) || badErr;
      double bw = h->GetBinWidth(1);
      printf("  %3.0f-%3.0f  | %4.0f |", lo, hi, h->Integral());
      if (bad) { printf("  COVARIANCE INVALID (status %d) -- bin dropped\n", covStat); continue; }
      TH_.push_back(0.5 * (lo + hi));
      for (int i = 0; i < NL; ++i) {
         double a = F->GetParameter(i) * SG[i] * std::sqrt(2 * TMath::Pi()) / bw;
         double e = F->GetParError(i) * SG[i] * std::sqrt(2 * TMath::Pi()) / bw;
         if (a > 0 && e < std::sqrt(a)) e = std::sqrt(a);          // floor at sqrt(N)
         Y[i].push_back(a); EY[i].push_back(e);
         printf(" %5.0f+-%3.0f |", a, e);
      }
      double aN = F->GetParameter(NL) * fa::sgN * std::sqrt(2 * TMath::Pi()) / bw;
      double aP = F->GetParameter(NL + 1) * (hps ? hps->Integral(hps->FindBin(exLo), hps->FindBin(exHi)) : 0);
      Y[NL].push_back(aN); EY[NL].push_back(std::sqrt(std::max(1.0, aN)));
      Y[NL + 1].push_back(aP); EY[NL + 1].push_back(std::sqrt(std::max(1.0, aP)));
      printf(" %5.0f | %6.0f | %6.2f\n", aN, aP,
             r.Get() && r->Ndf() > 0 ? r->Chi2() / r->Ndf() : -1);

      // ---- draw the bin with EVERY component ------------------------------------------------
      cf->cd(b + 1);
      gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
      h->GetXaxis()->SetLabelSize(0.05); h->GetYaxis()->SetLabelSize(0.05);
      h->GetXaxis()->SetTitleSize(0.055); h->GetYaxis()->SetTitleSize(0.055);
      h->GetXaxis()->SetTitleOffset(1.0); h->GetYaxis()->SetTitleOffset(1.1);
      h->SetTitle(Form("%.0f-%.0f deg", lo, hi));
      h->SetLineColor(kBlack); h->Draw("hist");
      F->SetLineColor(kRed); F->SetLineWidth(2); F->SetNpx(500); F->DrawCopy("same");
      for (int i = 0; i < NL; ++i) {          // one gaussian per level
         auto *g = new TF1(Form("g%d_%d", b, i), "gaus", exLo, exHi);
         g->SetParameters(F->GetParameter(i), E[i] + shift, SG[i]);
         g->SetLineColor(COL[i]); g->SetLineStyle(2); g->SetLineWidth(1); g->SetNpx(400); g->DrawCopy("same");
      }
      auto *gn = new TF1(Form("gn%d", b), "gaus", exLo, exHi);   // 14N (p,n)
      gn->SetParameters(F->GetParameter(NL), fa::muN, fa::sgN);
      gn->SetLineColor(COL[7]); gn->SetLineStyle(2); gn->SetLineWidth(2); gn->SetNpx(400); gn->DrawCopy("same");
      auto *gp = new TGraph(), *gb = new TGraph();               // continuum and background
      for (int k = 1; k <= h->GetNbinsX(); ++k) {
         double x = h->GetBinCenter(k);
         if (hps) gp->SetPoint(gp->GetN(), x, F->GetParameter(NL + 1) * hps->Interpolate(x));
         gb->SetPoint(gb->GetN(), x, F->GetParameter(NL + 2) + F->GetParameter(NL + 3) * x);
      }
      gp->SetLineColor(COL[8]); gp->SetLineWidth(2); gp->Draw("L same");
      gb->SetLineColor(COL[9]); gb->SetLineStyle(3); gb->SetLineWidth(2); gb->Draw("L same");
   }
   cf->SaveAs(here + "/plots/fit_angles_ps_fits_" + tag + ".png");

   // ---- acceptance and 1/sin(theta) -----------------------------------------------------------
   // Same convention as exc_angdist_C14.C: dsigma/dOmega ~ yield / acceptance(theta_cm) / sin(theta),
   // bins with acceptance <= 0.05 dropped. The multiplet uses the 6.094 acceptance and the 8.317
   // state its OWN (a dedicated simulation exists at that excitation energy); the 14N blend and the
   // continuum are NOT 14C levels, so no level acceptance applies and they are left raw.
   TH1D *accEx1 = nullptr, *accEx8 = nullptr;
   if (doAcc) {
      TFile *f1 = TFile::Open(accDir + "acceptance_merged_ex1.root");
      TFile *f8 = TFile::Open(accDir + "acceptance_merged_ex8.root");
      if (f1 && !f1->IsZombie()) { accEx1 = (TH1D *)f1->Get("hAcc_ex1_sum"); if (accEx1) accEx1->SetDirectory(nullptr); }
      if (f8 && !f8->IsZombie()) { accEx8 = (TH1D *)f8->Get("hAcc_ex8_sum"); if (accEx8) accEx8->SetDirectory(nullptr); }
      if (!accEx1) printf("\033[1;31m  no ex1 acceptance in %s -- yields left RAW\033[0m\n", accDir.Data());
   }
   std::vector<std::vector<double>> D(NL + 2), ED(NL + 2);
   std::vector<double> THd;
   if (accEx1) {
      printf("\n  acceptance-corrected, /sin(theta):\n  theta_cm |");
      for (int i = 0; i < NL; ++i) printf(" %10s |", JP[i]);
      printf("   acc\n");
      for (size_t k = 0; k < TH_.size(); ++k) {
         double c = TH_[k], sn = std::sin(c * TMath::DegToRad());
         double A1 = accEx1->GetBinContent(accEx1->FindBin(c));
         double A8 = accEx8 ? accEx8->GetBinContent(accEx8->FindBin(c)) : A1;
         if (A1 <= 0.05 || sn <= 1e-3) { printf("  %6.1f   acceptance %.3f -- DROPPED\n", c, A1); continue; }
         THd.push_back(c);
         printf("  %6.1f   |", c);
         for (int i = 0; i < NL + 2; ++i) {
            double A = (i == NL - 1 && accEx8) ? A8 : A1;      // 8.317 is the last fitted level
            if (i >= NL) A = 1.0;                              // 14N and continuum: no level acceptance
            D[i].push_back(Y[i][k] / A / sn);
            ED[i].push_back(EY[i][k] / A / sn);
            if (i < NL) printf(" %10.0f |", Y[i][k] / A / sn);
         }
         printf("  %.3f\n", A1);
      }
   }

   // ---- the angular distributions ------------------------------------------------------------
   auto *ca = new TCanvas("cang", "angular distributions", 1400, 900);
   ca->Divide(3, (NL + 4) / 3);
   std::vector<const char *> NM(NL + 2);
   for (int i = 0; i < NL; ++i) NM[i] = JP[i];
   NM[NL] = "14N (p,n) blend"; NM[NL + 1] = "1n continuum";
   for (int i = 0; i < NL + 2; ++i) {
      ca->cd(i + 1);
      bool corr = accEx1 && !THd.empty() && D[i].size() == THd.size();
      auto *g = corr ? new TGraphErrors(THd.size(), &THd[0], &D[i][0], nullptr, &ED[i][0])
                     : new TGraphErrors(TH_.size(), &TH_[0], &Y[i][0], nullptr, &EY[i][0]);
      g->SetTitle(Form("%s;#theta_{cm} [deg];%s", NM[i],
                       corr ? "yield / acc / sin#theta  [arb.]" : "counts per bin"));
      g->SetMarkerStyle(20); g->SetMarkerColor(COL[i < 7 ? i : 7 + (i - NL)]); g->SetLineColor(COL[i < 7 ? i : 7 + (i - NL)]);
      // LOG y. Points at or below zero cannot be drawn on a log axis, and silently vanishing is
      // worse than being told: they are removed from the graph and counted in the title.
      int nzero = 0;
      for (int k = g->GetN() - 1; k >= 0; --k)
         if (g->GetY()[k] <= 0) { g->RemovePoint(k); ++nzero; }
      double ymax = 0, ymin = 1e30;
      for (int k = 0; k < g->GetN(); ++k) {
         ymax = std::max(ymax, g->GetY()[k] + g->GetEY()[k]);
         ymin = std::min(ymin, g->GetY()[k]);
      }
      if (g->GetN() == 0) { continue; }
      if (nzero) g->SetTitle(Form("%s  [%d bin%s at zero not shown];#theta_{cm} [deg];%s", NM[i],
                                  nzero, nzero > 1 ? "s" : "",
                                  corr ? "yield / acc / sin#theta  [arb.]" : "counts per bin"));
      gPad->SetLogy();
      g->SetMinimum(0.5 * std::max(1e-3, ymin));
      g->SetMaximum(2.0 * ymax);
      g->Draw("AP");
   }
   ca->SaveAs(here + "/plots/fit_angles_ps_dist_" + tag + ".png");
   printf("\n  wrote plots/fit_angles_ps_fits_%s.png and plots/fit_angles_ps_dist_%s.png\n", tag.Data(), tag.Data());
}
