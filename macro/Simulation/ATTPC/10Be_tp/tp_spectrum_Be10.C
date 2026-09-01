/// @file tp_spectrum_Be10.C
/// @brief The 12Be excitation-energy spectrum the AT-TPC would actually measure: the four bound
///        levels SUMMED with their relative populations, and the test of whether the 0+_2 is
///        recoverable from it.
///
///   root -b -q 'tp_spectrum_Be10.C("/mnt/f/Be10_tp")'
///   root -b -q 'tp_spectrum_Be10.C("/mnt/f/Be10_tp","b700_2mm",0.2,2.0,45.0)'
///
/// POPULATIONS. 0+ g.s., 2+ 2.109 and 1- 2.715 enter with equal weight; the intruder 0+_2 at
/// 2.251 enters 5x weaker (w2251 = 0.2), as asked. The suppression is applied HERE and not in the
/// generation: every level was simulated at full statistics so that its acceptance and its
/// resolution are each measured to the same precision, and the population is a multiplicative
/// weight on a spectrum, not a property of the detector. Generating a fifth of the events would
/// have made the hardest level the worst-measured one.
///
/// The weight actually applied to level i is  pop_i / nGen_i , where nGen_i is the number of
/// reactions GENERATED for that level (read from the acceptance file, not assumed equal). That way
/// each level's own acceptance survives into the sum: the histogram is proportional to what a
/// detector would count, not to what was thrown.
///
/// THE TEST. A resolution number does not by itself say whether a 142 keV doublet is separable
/// when one member is five times weaker. So this macro does what an experiment would have to do:
///   FIT A  four gaussians, everything free from a sensible start -- does the fit find the 0+_2,
///          and does it return its input area ratio of 0.2?
///   FIT B  three gaussians, the doublet replaced by ONE peak -- the null hypothesis.
/// Both are Poisson-likelihood fits (option "L"), and they are compared on the LIKELIHOOD, not on
/// chi2, because chi2 is not the quantity being minimised. A large drop from B to A means the
/// spectrum demands the extra state; a small one means the 0+_2 is invisible at this resolution
/// however good the fit to it happens to look.

#include <algorithm>
#include <map>
#include <vector>

static const int NL = 4;
static const char *LVL[NL] = {"gs", "ex2109", "ex2251", "ex2715"};
static const double LEX[NL] = {0.0, 2.109, 2.251, 2.715};
static const char *LJP[NL] = {"0+", "2+", "0+_2", "1-"};

static TString sp_find(const TString &dir, const TString &cfg, const char *lvl, const char *pre)
{
   TString f = gSystem->GetFromPipe(
      TString::Format("ls %s/%s/%s_%s_s*_%s.root 2>/dev/null | head -1", dir.Data(), cfg.Data(), pre, lvl, cfg.Data()));
   return f.Strip(TString::kBoth);
}

/// number of reactions GENERATED for this level, from the acceptance file's truth denominator
static double sp_ngen(const TString &dir, const TString &cfg, const char *lvl)
{
   TString fa = sp_find(dir, cfg, lvl, "acceptance");
   if (fa.IsNull()) return -1;
   TFile *f = TFile::Open(fa);
   if (!f || f->IsZombie()) return -1;
   double n = -1;
   TIter nx(f->GetListOfKeys());
   while (auto *k = (TKey *)nx()) {
      TString nm = k->GetName();
      if (nm.BeginsWith("hGen_")) { TH1D *h = (TH1D *)f->Get(nm); if (h) n = h->Integral(); }
   }
   f->Close();
   return n;
}

void tp_spectrum_Be10(TString root = "/mnt/f/Be10_tp", TString cfg = "b285_attpc", Double_t w2251 = 0.2,
                      Double_t cmLo = 2.0, Double_t cmHi = 180.0, Double_t exLo = -2.0, Double_t exHi = 5.0,
                      Int_t nBins = 140, TString outDir = "")
{
   gStyle->SetOptStat(0);
   if (outDir.IsNull()) outDir = gSystem->pwd();
   const double pop[NL] = {1.0, 1.0, w2251, 1.0};

   TH1D *hSum = new TH1D("hSum", ";E_{x}(^{12}Be) [MeV];counts / bin", nBins, exLo, exHi);
   TH1D *hLev[NL];
   double nGen[NL], nUse[NL], wgt[NL];
   bool have = true;
   printf("\n=========== 10Be(t,p)12Be spectrum : %s / %s ===========\n", root.Data(), cfg.Data());
   printf("theta_cm %.0f-%.0f deg, populations 1 : 1 : %.2f : 1 (the 0+_2 suppressed %.1fx)\n", cmLo, cmHi, w2251,
          w2251 > 0 ? 1.0 / w2251 : 0.0);
   for (int l = 0; l < NL; ++l) {
      hLev[l] = new TH1D(Form("hL%d", l), "", nBins, exLo, exHi);
      hLev[l]->SetDirectory(nullptr);
      nGen[l] = sp_ngen(root, cfg, LVL[l]);
      nUse[l] = 0;
      wgt[l] = 0;
      TString f = sp_find(root, cfg, LVL[l], "exres");
      if (f.IsNull() || nGen[l] <= 0) {
         printf("\033[1;31mMISSING level %s in %s (exres '%s', nGen %.0f)\033[0m\n", LVL[l], cfg.Data(), f.Data(),
                nGen[l]);
         have = false;
         continue;
      }
      TFile *fr = TFile::Open(f);
      TTree *t = (TTree *)fr->Get("res");
      double exReco, cmTrue;
      t->SetBranchAddress("exReco", &exReco);
      t->SetBranchAddress("cmTrue", &cmTrue);
      wgt[l] = pop[l] / nGen[l];
      for (Long64_t i = 0; i < t->GetEntries(); ++i) {
         t->GetEntry(i);
         if (cmTrue < cmLo || cmTrue >= cmHi) continue;
         hLev[l]->Fill(exReco, wgt[l]);
         nUse[l] += 1;
      }
      fr->Close();
      hSum->Add(hLev[l]);
   }
   if (!have) { printf("incomplete configuration -- nothing to fit\n"); return; }

   printf("\n  %-8s %-6s %8s %8s %10s %12s\n", "level", "J^pi", "nGen", "nFit", "weight", "yield (rel)");
   double totYield = 0;
   for (int l = 0; l < NL; ++l) totYield += hLev[l]->Integral();
   for (int l = 0; l < NL; ++l)
      printf("  %-8s %-6s %8.0f %8.0f %10.3e %11.4f  (acceptance %.3f)\n", LVL[l], LJP[l], nGen[l], nUse[l], wgt[l],
             hLev[l]->Integral() / totYield, nUse[l] / nGen[l]);

   // Scale to a realistic number of counts so the fit sees Poisson statistics rather than weights.
   // 20000 total counts is the size of a decent (t,p) run and keeps the 0+_2 at ~1500 counts, so a
   // failure to find it is a RESOLUTION failure and not a statistics one -- which is the point.
   const double NTOT = 20000.0;
   if (hSum->Integral() > 0) hSum->Scale(NTOT / hSum->Integral());
   for (int l = 0; l < NL; ++l) hLev[l]->Scale(NTOT / totYield);
   // Poissonise, so the likelihood fit is fitting counts
   TH1D *hObs = (TH1D *)hSum->Clone("hObs");
   hObs->SetDirectory(nullptr);
   for (int b = 1; b <= hObs->GetNbinsX(); ++b) hObs->SetBinContent(b, gRandom->Poisson(hSum->GetBinContent(b)));
   hObs->SetTitle(Form("10Be(t,p)12Be, %s, #theta_{cm} %.0f-%.0f deg", cfg.Data(), cmLo, cmHi));

   // ---- FIT A : four gaussians ------------------------------------------------------------
   TF1 *fA = new TF1("fA", "gaus(0)+gaus(3)+gaus(6)+gaus(9)", exLo, exHi);
   double sig0 = 0.4;
   for (int l = 0; l < NL; ++l) {
      fA->SetParameter(3 * l + 0, hObs->GetMaximum() * (l == 2 ? 0.2 : 0.8));
      fA->SetParameter(3 * l + 1, LEX[l]);
      fA->SetParameter(3 * l + 2, sig0);
      fA->SetParLimits(3 * l + 0, 0, 1e7);
      fA->SetParLimits(3 * l + 1, LEX[l] - 1.0, LEX[l] + 1.0);
      fA->SetParLimits(3 * l + 2, 0.02, 3.0);
      fA->SetParName(3 * l + 0, Form("A_%s", LVL[l]));
      fA->SetParName(3 * l + 1, Form("mu_%s", LVL[l]));
      fA->SetParName(3 * l + 2, Form("sig_%s", LVL[l]));
   }
   auto rA = hObs->Fit(fA, "LSQ0", "", exLo, exHi);

   // ---- FIT B : three gaussians, the 2.109/2.251 doublet as ONE peak ----------------------
   TF1 *fB = new TF1("fB", "gaus(0)+gaus(3)+gaus(6)", exLo, exHi);
   const double bex[3] = {0.0, 2.15, 2.715};
   for (int l = 0; l < 3; ++l) {
      fB->SetParameter(3 * l + 0, hObs->GetMaximum() * 0.8);
      fB->SetParameter(3 * l + 1, bex[l]);
      fB->SetParameter(3 * l + 2, sig0);
      fB->SetParLimits(3 * l + 0, 0, 1e7);
      fB->SetParLimits(3 * l + 1, bex[l] - 1.0, bex[l] + 1.0);
      fB->SetParLimits(3 * l + 2, 0.02, 3.0);
   }
   auto rB = hObs->Fit(fB, "LSQ0+", "", exLo, exHi);

   printf("\n---- FIT A : four gaussians (every parameter listed, and flagged if it sits on a bound) ----\n");
   printf("  %-12s %12s %12s %10s\n", "parameter", "value", "error", "");
   for (int p = 0; p < fA->GetNpar(); ++p) {
      double lo, hi;
      fA->GetParLimits(p, lo, hi);
      double v = fA->GetParameter(p);
      bool rail = (hi > lo) && (std::fabs(v - lo) < 1e-6 * std::max(1.0, std::fabs(lo)) ||
                                std::fabs(v - hi) < 1e-6 * std::max(1.0, std::fabs(hi)));
      printf("  %-12s %12.5f %12.5f %10s\n", fA->GetParName(p), v, fA->GetParError(p), rail ? "<-- RAILED" : "");
   }
   // areas, from A*sigma*sqrt(2pi) / bin width
   double bw = hObs->GetBinWidth(1);
   double areaA[NL];
   for (int l = 0; l < NL; ++l)
      areaA[l] = fA->GetParameter(3 * l) * fA->GetParameter(3 * l + 2) * std::sqrt(2 * TMath::Pi()) / bw;
   printf("\n  %-8s %10s %10s %12s %12s %12s\n", "level", "Ex in", "Ex fit", "sigma fit", "area", "area/area(2+)");
   for (int l = 0; l < NL; ++l)
      printf("  %-8s %10.3f %10.3f %12.3f %12.1f %12.3f\n", LVL[l], LEX[l], fA->GetParameter(3 * l + 1),
             fA->GetParameter(3 * l + 2), areaA[l], areaA[1] > 0 ? areaA[l] / areaA[1] : NAN);
   printf("  the 0+_2 was PUT IN at area/area(2+) = %.3f; the fit returns %.3f\n", w2251,
          areaA[1] > 0 ? areaA[2] / areaA[1] : NAN);

   double lA = rA.Get() ? rA->MinFcnValue() : NAN;
   double lB = rB.Get() ? rB->MinFcnValue() : NAN;
   printf("\n---- FIT A vs FIT B : does the spectrum DEMAND the 0+_2? ----\n");
   printf("  compared on the LIKELIHOOD (option \"L\" minimises -2lnL; chi2 is not the objective)\n");
   printf("  A (4 peaks, %d par) : -2lnL = %.2f\n", fA->GetNpar(), lA);
   printf("  B (3 peaks, %d par) : -2lnL = %.2f\n", fB->GetNpar(), lB);
   double dl = lB - lA;
   printf("  \033[1;3%dm Delta(-2lnL) = %.2f for 3 extra parameters  -> %s\033[0m\n", dl > 11.3 ? 2 : 1, dl,
          dl > 11.3 ? "the extra state is DEMANDED (>3 sigma)"
                    : (dl > 6.3 ? "marginal (2-3 sigma)" : "the 0+_2 is NOT recoverable from this spectrum"));
   printf("  (11.3 and 6.3 are the 3- and 2-sigma points of a chi2 with 3 degrees of freedom)\n");

   // ---- figure -----------------------------------------------------------------------------
   TCanvas *c = new TCanvas("cSpec", "spectrum", 1200, 800);
   c->Divide(1, 2);
   c->cd(1);
   gPad->SetPad(0, 0.35, 1, 1);
   hObs->SetLineColor(kBlack);
   hObs->SetMarkerStyle(20);
   hObs->SetMarkerSize(0.5);
   hObs->GetXaxis()->SetRangeUser(exLo, exHi);
   hObs->Draw("E");
   fA->SetLineColor(kRed);
   fA->SetNpx(1000);
   fA->Draw("same");
   int col[NL] = {kBlue, kGreen + 2, kMagenta, kOrange + 7};
   TF1 *comp[NL];
   for (int l = 0; l < NL; ++l) {
      comp[l] = new TF1(Form("c%d", l), "gaus", exLo, exHi);
      comp[l]->SetParameters(fA->GetParameter(3 * l), fA->GetParameter(3 * l + 1), fA->GetParameter(3 * l + 2));
      comp[l]->SetLineColor(col[l]);
      comp[l]->SetLineStyle(2);
      comp[l]->SetNpx(1000);
      comp[l]->Draw("same");
   }
   TLegend *leg = new TLegend(0.62, 0.55, 0.98, 0.90);
   leg->AddEntry(hObs, Form("sum, %s", cfg.Data()), "lep");
   leg->AddEntry(fA, "4-gaussian fit", "l");
   for (int l = 0; l < NL; ++l)
      leg->AddEntry(comp[l], Form("%s  %.3f MeV  (w %.2f)", LJP[l], LEX[l], pop[l]), "l");
   leg->Draw();
   c->cd(2);
   gPad->SetPad(0, 0, 1, 0.35);
   TH1D *hR = (TH1D *)hObs->Clone("hR");
   hR->SetDirectory(nullptr);
   for (int b = 1; b <= hR->GetNbinsX(); ++b) {
      double m = fA->Eval(hR->GetBinCenter(b));
      double e = std::sqrt(std::max(1.0, m));
      hR->SetBinContent(b, (hObs->GetBinContent(b) - m) / e);
      hR->SetBinError(b, 1);
   }
   hR->SetTitle(";E_{x} [MeV];(data - fit A)/#sigma");
   hR->GetYaxis()->SetRangeUser(-5, 5);
   hR->Draw("E");
   TLine *l0 = new TLine(exLo, 0, exHi, 0);
   l0->SetLineColor(kRed);
   l0->Draw();
   TString png = TString::Format("%s/spectrum_Be10tp_%s_cm%03.0f-%03.0f.png", outDir.Data(), cfg.Data(), cmLo, cmHi);
   c->SaveAs(png);
   printf("\nwrote %s\n", png.Data());
   printf("spectrum done\n\n");
}
