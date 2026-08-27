/// @file spec_states_C14.C
/// @brief GLOBAL fit of the Ex spectrum over one theta_cm window, with different numbers of states.
///
/// One spectrum, summed over theta_cm 20-60 deg, fitted repeatedly with a growing level set. The
/// question is how many states the SPECTRUM itself demands -- separate from what the angular
/// distributions do with them.
///
/// Peak positions are fixed at literature + the global shift and widths at sig0 + dSig*(E-6.094),
/// exactly as fit_states_ps_C14.C does; only the amplitudes, the 14N blend, the 1n continuum and a
/// linear background float. Likelihood fit; the Pearson chi2/ndf is reported alongside.
///
///   root -b -q 'spec_states_C14.C()'
namespace ss {
int NL = 5;
double E[8], SG[8];
double shift = -0.011, muN = 9.178, sgN = 0.296;
TH1D *hps = nullptr;
double model(double *x, double *p)
{
   double s = p[NL + 2] + p[NL + 3] * x[0];
   for (int i = 0; i < NL; ++i)
      s += p[i] * std::exp(-0.5 * std::pow((x[0] - (E[i] + shift)) / SG[i], 2));
   s += p[NL] * std::exp(-0.5 * std::pow((x[0] - muN) / sgN, 2));
   if (hps) s += p[NL + 1] * hps->Interpolate(x[0]);
   return s;
}
} // namespace ss

void spec_states_C14(TString cache = "plots/proton_kin_cat5_tc.root",
                     Double_t cmLo = 20, Double_t cmHi = 60,
                     Double_t exLo = 5.0, Double_t exHi = 10.4, Double_t shift = -0.011,
                     Double_t sig0 = 0.132, Double_t dSig = 0.0123,
                     Double_t vzLo = 10, Double_t vzHi = 490, Double_t chi2Cut = 5.0)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TFile *fd = TFile::Open(here + "/" + cache);
   TNtuple *t = fd && !fd->IsZombie() ? (TNtuple *)fd->Get("pk") : nullptr;
   TFile *fp = TFile::Open(here + "/plots/phasespace_1n_C14.root");
   if (!t) { printf("\033[1;31mno cache\033[0m\n"); return; }
   if (fp && !fp->IsZombie()) { ss::hps = (TH1D *)fp->Get("hps"); if (ss::hps) ss::hps->SetDirectory(nullptr); }
   ss::shift = shift;

   const int nb = (int)std::lround((exHi - exLo) / 0.05);
   auto *h = new TH1D("hspec", "", nb, exLo, exHi);
   h->Sumw2();
   t->Draw("ex>>hspec",
           Form("chi2ndf<%g && thcm>=%g && thcm<%g && vertexz>=%g && vertexz<=%g",
                chi2Cut, cmLo, cmHi, vzLo, vzHi), "goff");
   printf("\n  Ex spectrum, theta_cm %.0f-%.0f deg : %.0f counts in %.1f-%.1f MeV\n",
          cmLo, cmHi, h->Integral(), exLo, exHi);

   // growing level sets. Each row adds ONE state to the row above.
   const int NS = 5;
   const char *setName[NS] = {"4 states", "5 states  (+8.317)", "6 states  (+6.589)",
                              "7 states  (+6.903)", "8 states  (+9.80)"};
   std::vector<std::vector<double>> ES = {
      {6.091, 6.728, 7.012, 7.341},
      {6.091, 6.728, 7.012, 7.341, 8.317},
      {6.091, 6.589, 6.728, 7.012, 7.341, 8.317},
      {6.091, 6.589, 6.728, 6.903, 7.012, 7.341, 8.317},
      {6.091, 6.589, 6.728, 6.903, 7.012, 7.341, 8.317, 9.801}};

   auto *cv = new TCanvas("cs", "", 1600, 1000); cv->Divide(2, 3);
   printf("\n  %-22s %5s %9s %9s | %10s %10s\n", "level set", "npar", "chi2/ndf", "-lnL", "6.728 area", "7.012 area");
   printf("  %s\n", TString('-', 76).Data());
   double prevC2 = 0;
   for (int s = 0; s < NS; ++s) {
      ss::NL = (int)ES[s].size();
      for (int i = 0; i < ss::NL; ++i) { ss::E[i] = ES[s][i]; ss::SG[i] = sig0 + dSig * (ES[s][i] - 6.094); }
      const int NPAR = ss::NL + 4;
      auto *F = new TF1(Form("F%d", s), ss::model, exLo, exHi, NPAR);
      for (int i = 0; i < ss::NL; ++i) {
         F->SetParameter(i, std::max(1.0, h->GetBinContent(h->FindBin(ES[s][i] + shift)) * 0.5));
         F->SetParLimits(i, 0, 1e5);
      }
      F->SetParameter(ss::NL, 20);      F->SetParLimits(ss::NL, 0, 1e5);
      F->SetParameter(ss::NL + 1, h->Integral() * 0.05); F->SetParLimits(ss::NL + 1, 0, 1e7);
      F->SetParameter(ss::NL + 2, 2);   F->SetParameter(ss::NL + 3, 0);
      TFitResultPtr r = h->Fit(F, "RQNSL");
      for (int pass = 0; pass < 3 && (int)r->Status() != 0; ++pass) r = h->Fit(F, "RQNSL");
      double c2 = 0; int nfree = 0;
      for (int b = 1; b <= nb; ++b) {
         double m = F->Eval(h->GetBinCenter(b)), o = h->GetBinContent(b);
         if (m > 0) { c2 += std::pow(o - m, 2) / m; ++nfree; }
      }
      int ndf = nfree - NPAR;
      double bw = h->GetBinWidth(1), sq = std::sqrt(2 * TMath::Pi());
      auto area = [&](double Ei) {
         for (int i = 0; i < ss::NL; ++i)
            if (std::fabs(ES[s][i] - Ei) < 0.01) return F->GetParameter(i) * ss::SG[i] * sq / bw;
         return 0.0; };
      printf("  %-22s %5d %9.2f %9.1f | %10.0f %10.0f\n", setName[s], NPAR,
             ndf > 0 ? c2 / ndf : 0, r->MinFcnValue(), area(6.728), area(7.012));
      if (s) { double d = prevC2 - c2 / ndf; (void)d; }
      prevC2 = c2 / ndf;

      cv->cd(s + 1); gPad->SetGridx(); gPad->SetGridy();
      auto *hc = (TH1D *)h->Clone(Form("hc%d", s));
      hc->SetTitle(Form("%s   #chi^{2}/ndf = %.2f;E_{x} [MeV];counts / 50 keV",
                        setName[s], ndf > 0 ? c2 / ndf : 0));
      hc->SetMarkerStyle(20); hc->SetMarkerSize(0.5); hc->SetLineColor(kBlack);
      hc->Draw("E");
      // A TF1 built on a C++ function evaluates LAZILY: it re-reads ss::NL / ss::E / ss::SG when the
      // canvas is finally drawn, by which time the loop has moved on and every panel renders with
      // the LAST state set. Snapshot the curve into a graph now, while the globals are still this
      // panel's.
      auto *gt = new TGraph();
      for (double x = exLo; x <= exHi; x += 0.01) gt->SetPoint(gt->GetN(), x, F->Eval(x));
      gt->SetLineColor(kRed + 1); gt->SetLineWidth(3); gt->Draw("L same");
      // components
      const int cc[9] = {kAzure + 2, kGreen + 3, kMagenta + 2, kOrange + 8, kCyan + 2,
                         kViolet + 1, kSpring - 6, kPink + 7, kGray + 2};
      for (int i = 0; i < ss::NL; ++i) {
         auto *g = new TF1(Form("g%d_%d", s, i), "gaus", exLo, exHi);
         g->SetParameters(F->GetParameter(i), ES[s][i] + shift, ss::SG[i]);
         g->SetLineColor(cc[i]); g->SetLineStyle(2); g->SetLineWidth(2); g->SetNpx(400);
         g->DrawCopy("same");
      }
      TLatex tx; tx.SetNDC(); tx.SetTextSize(0.040);
      tx.DrawLatex(0.60, 0.84, Form("6.728 : %.0f", area(6.728)));
      tx.DrawLatex(0.60, 0.79, Form("7.012 : %.0f", area(7.012)));
   }
   cv->cd(6);
   auto *pt = new TPaveText(0.03, 0.05, 0.97, 0.95); pt->SetBorderSize(0); pt->SetFillStyle(0);
   pt->SetTextAlign(12); pt->SetTextSize(0.048);
   pt->AddText(Form("#bf{E_{x} spectrum, #theta_{cm} = %.0f-%.0f deg}", cmLo, cmHi));
   pt->AddText("");
   pt->AddText("each panel adds ONE state to the panel before it");
   pt->AddText("positions fixed at literature + shift,");
   pt->AddText("widths at #sigma_{0} + d#sigma(E - 6.094)");
   pt->AddText("");
   pt->AddText("dashed = individual levels, red = total");
   pt->Draw();
   TString out = "/home/yassid/a1954_analysis_runs/2026-08-25_C14_catima_refit/plots/06_ptolemy/";
   cv->SaveAs(out + "26_spectrum_state_sets.png");
   printf("\n  wrote %s26_spectrum_state_sets.png\n\n", out.Data());
}
