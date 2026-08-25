/// @file fit_states_ps_C14.C
/// @brief ONE fit of the 14C(p,p') spectrum from 6 to 10.4 MeV, levels PLUS the 1n-emission
///        continuum. Sn(14C) = 8.176 MeV sits inside the window, so a fit with no continuum term
///        must absorb it into the peak areas or a straight line.
///
/// Recipe from the 16C(p,d) phase-space fit (fit_pd_ps.C), every point of which was paid for:
///   * GAUSSIANS, never Breit-Wigners -- a BW collapses to a delta spike or runs to its bound and
///     acts as extra background.
///   * WIDTHS FIXED at the resolution model sigma(Ex) = 0.132 + 0.0123 (Ex - 6.094). Widths are
///     resolution and cannot be allowed to depend on the background model.
///   * POSITIONS FROZEN at literature + the measured shift. Given a window they slide and trade.
///   * POISSON LIKELIHOOD ("L"), not chi2: chi2 discards empty bins, so above 9.5 where bins are
///     sparse the continuum would be unconstrained and invent counts, and Neyman chi2 biases
///     amplitudes low as bins empty.
///   * Errors floored at sqrt(N): Minuit returns absurdly small errors on small areas.
/// The 14N (p,n) structure at ~9.18 is a DIFFERENT REACTION, not a 14C level and not continuum;
/// it gets its own gaussian at its measured width (0.296, twice the resolution -- a blend).
///
///   root -b -q 'fit_states_ps_C14.C("plots/proton_kin_cat5_tc.root")'
void fit_states_ps_C14(TString cache = "plots/proton_kin_cat5_tc.root", Double_t chi2Cut = 5.0,
                       Double_t exLo = 6.0, Double_t exHi = 10.4, Double_t shift = -0.011,
                       Double_t sig0 = 0.132, Double_t dSig = 0.0123, Double_t muN = 9.178,
                       Double_t sgN = 0.296, Int_t nb = 88, Bool_t usePS = kTRUE, TString tag = "cat5",
                       // ANCHOR THE BACKGROUND IN THE GAP INSTEAD OF FLOATING IT.
                       // There are no 14C levels between the g.s. and 6.094, so 3.0-5.0 MeV is a
                       // DIRECT measurement of the background: it falls off the g.s. tail (12.1
                       // counts/50 keV at 1-2 MeV) and flattens at ~3.2 by 3 MeV. A straight line
                       // fitted only above 6 MeV cannot know that and came out at ~13/bin at 6 MeV,
                       // four times the measured floor -- absorbing peak area and making the yields
                       // depend on the background model. Fixing it from the gap breaks that
                       // degeneracy, which is the same trick the 16C(p,d) phase-space fit uses.
                       // MEASURED AND SET ASIDE (2026-08-26). The gap level is 3.17 counts per
                       // 50 keV bin, and at 9.6-10.4 the data agrees with it (4.31). But at 7.7-8.1
                       // -- BELOW Sn, so not continuum -- the data sits at 7.62, 2.4x the gap, and
                       // no gaussian reaches there (7.9 is 3.7 sigma from the 7.341 peak). That is
                       // the HIGH-SIDE TAIL of the peaks: the g.s. shows the same thing, 12.05 per
                       // bin at 1-2 MeV against the same 3.17 floor, which is 7-14 sigma out on a
                       // sigma = 0.144 peak. Pinning the background at the gap therefore makes the
                       // fit WORSE (chi2/ndf 1.50 -> 1.94) and collapses the continuum, because the
                       // free line had been standing in for the tail. The real fix is a lineshape
                       // with a tail, not a better background. Until then the background is left
                       // FREE -- it is 9 % of the window and Yassid's call is that it is small
                       // enough to set aside. Pass bgGapLo > 0 to reproduce the pinned version.
                       Double_t bgGapLo = -1, Double_t bgGapHi = 5.0)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TFile *fd = TFile::Open(here + "/" + cache);
   TNtuple *t = (TNtuple *)fd->Get("pk");
   auto *h = new TH1D("hPSfit", "^{14}C(p,p') 6-10.4 MeV, levels + 1n continuum;E_{x} [MeV];counts",
                      nb, exLo, exHi);
   h->Sumw2();
   t->Draw("ex>>hPSfit", TString::Format("chi2ndf<%g", chi2Cut), "goff");

   static TH1D *hps = nullptr;
   if (usePS) {
      TFile *fp = TFile::Open(here + "/plots/phasespace_1n_C14.root");
      if (!fp || fp->IsZombie()) { printf("\033[1;31mrun phasespace_1n_C14.C first\033[0m\n"); return; }
      hps = (TH1D *)fp->Get("hPS");
      hps->SetDirectory(nullptr);
      hps->Scale(1.0 / hps->Integral());     // shape only; one free scale carries the normalisation
   }
   // 14C levels: literature energies, shifted by the measured common shift
   const int NL = 7;
   double E[NL] = {6.091, 6.589, 6.728, 6.903, 7.012, 7.341, 8.317};
   const char *JP[NL] = {"1-", "0+_2", "3-_1", "0-_1", "2+_1", "2-_1", "2+ (8.317)"};
   double SG[NL];
   for (int i = 0; i < NL; ++i) SG[i] = sig0 + dSig * (E[i] - 6.094);

   const int NPAR = NL + 1 /*14N*/ + 1 /*PS*/ + 2 /*bg*/;
   auto model = [&](double *x, double *p) {
      double s = p[NL + 2] + p[NL + 3] * x[0];
      for (int i = 0; i < NL; ++i) s += p[i] * std::exp(-0.5 * std::pow((x[0] - (E[i] + shift)) / SG[i], 2));
      s += p[NL] * std::exp(-0.5 * std::pow((x[0] - muN) / sgN, 2));
      if (usePS && hps) s += p[NL + 1] * hps->Interpolate(x[0]);
      return s;
   };
   TF1 F("F", model, exLo, exHi, NPAR);
   for (int i = 0; i < NL; ++i) { F.SetParameter(i, std::max(1.0, h->GetBinContent(h->FindBin(E[i] + shift)) * 0.5));
                                  F.SetParLimits(i, 0, 1e5); F.SetParName(i, JP[i]); }
   F.SetParameter(NL, 50); F.SetParLimits(NL, 0, 1e5); F.SetParName(NL, "14N(p,n)");
   F.SetParameter(NL + 1, usePS ? h->Integral() * 0.3 : 0);
   if (usePS) F.SetParLimits(NL + 1, 0, 1e7); else F.FixParameter(NL + 1, 0);  // limits BEFORE fix elsewhere
   F.SetParName(NL + 1, "1n phase space");
   F.SetParameter(NL + 2, 1); F.SetParameter(NL + 3, 0);
   F.SetParName(NL + 2, "bg0"); F.SetParName(NL + 3, "bg1");
   double bgLevel = -1;
   if (bgGapLo > 0) {
      // measure it on the same cache, same bin width as the fit histogram
      auto *hg = new TH1D("hGap", "", (int)((bgGapHi - bgGapLo) / h->GetBinWidth(1) + 0.5), bgGapLo, bgGapHi);
      t->Draw("ex>>hGap", TString::Format("chi2ndf<%g", chi2Cut), "goff");
      bgLevel = hg->Integral() / hg->GetNbinsX();
      // NOTE: no SetParLimits after this. SetParLimits AFTER FixParameter silently RELEASES the
      // parameter, and the tell is that the output is byte-identical to the unfixed fit.
      F.FixParameter(NL + 2, bgLevel);
      F.FixParameter(NL + 3, 0.0);
      printf("\n  background FIXED from the gap %.1f-%.1f MeV: %.2f counts per %.0f keV bin"
             "  (%.0f counts in the gap)\n", bgGapLo, bgGapHi, bgLevel, 1000 * h->GetBinWidth(1),
             hg->Integral());
   }
   TFitResultPtr r = h->Fit(&F, "RQNSL");         // L = Poisson likelihood

   printf("\n\033[1;33m===== 6-10.4 MeV, %s 1n continuum =====\033[0m\n", usePS ? "WITH" : "WITHOUT");
   double bw = h->GetBinWidth(1), tot = 0;
   for (int i = 0; i < NL; ++i) {
      double a = F.GetParameter(i) * SG[i] * std::sqrt(2 * TMath::Pi()) / bw;
      double e = F.GetParError(i) * SG[i] * std::sqrt(2 * TMath::Pi()) / bw;
      if (a > 0 && e < std::sqrt(a)) e = std::sqrt(a);        // floor the error at sqrt(N)
      printf("   %-12s %6.3f   area %7.0f +- %5.0f%s\n", JP[i], E[i] + shift, a, e,
             F.GetParameter(i) < 1e-6 ? "   <-- COLLAPSED" : "");
      tot += a;
   }
   double aN = F.GetParameter(NL) * sgN * std::sqrt(2 * TMath::Pi()) / bw;
   printf("   %-12s %6.3f   area %7.0f            <-- NOT 14C, (p,n) to 14N\n", "14N blend", muN, aN);
   if (usePS) {
      double aps = F.GetParameter(NL + 1) * (hps ? hps->Integral(hps->FindBin(exLo), hps->FindBin(exHi)) : 0);
      printf("   %-12s          area %7.0f  (%.0f %% of all counts in the window)\n", "1n continuum", aps,
             100.0 * aps / std::max(1.0, h->Integral()));
      printf("        of which above S_n = 8.176 : %.1f %%\n",
             100.0 * hps->Integral(hps->FindBin(8.176), hps->FindBin(exHi))
                   / std::max(1e-9, hps->Integral(hps->FindBin(exLo), hps->FindBin(exHi))));
   }
   {
      double bgI = 0;
      for (int b = 1; b <= h->GetNbinsX(); ++b)
         bgI += F.GetParameter(NL + 2) + F.GetParameter(NL + 3) * h->GetBinCenter(b);
      printf("   %-12s          area %7.0f  (%.0f %% of all counts in the window)\n", "flat bg", bgI,
             100.0 * bgI / std::max(1.0, h->Integral()));
   }
   printf("   14C peak total %.0f    counts in window %.0f\n", tot, h->Integral());
   printf("   chi2/ndf = %.2f / %d = %.2f   (Pearson, reported alongside the likelihood fit)\n",
          r.Get() ? r->Chi2() : -1, r.Get() ? r->Ndf() : 0,
          r.Get() && r->Ndf() > 0 ? r->Chi2() / r->Ndf() : -1);

   auto *c = new TCanvas("cps", "", 1000, 700);
   h->SetLineColor(kBlack); h->Draw("hist");
   F.SetLineColor(kRed + 1); F.SetNpx(600); F.DrawCopy("same");
   // Draw the continuum and the flat background SEPARATELY. Drawing their SUM invites exactly the
   // right question -- "why is there phase space below Sn?" -- when the answer is that there is
   // none: the 1n term is zero below 7.6 and carries 0.9 % of its weight below Sn = 8.176, all of
   // that being resolution smearing across the threshold.
   auto *gb = new TGraph();
   for (int b = 1; b <= h->GetNbinsX(); ++b) {
      double x = h->GetBinCenter(b);
      gb->SetPoint(gb->GetN(), x, F.GetParameter(NL + 2) + F.GetParameter(NL + 3) * x);
   }
   gb->SetLineColor(kGreen + 2); gb->SetLineWidth(2); gb->SetLineStyle(3); gb->Draw("L same");
   TLegend *lg = new TLegend(0.58, 0.68, 0.89, 0.88);
   lg->SetBorderSize(0); lg->SetFillStyle(0);
   lg->AddEntry(h, "data", "l");
   lg->AddEntry(&F, "total fit", "l");
   if (usePS && hps) {
      auto *g = new TGraph();
      for (int b = 1; b <= h->GetNbinsX(); ++b) {
         double x = h->GetBinCenter(b);
         g->SetPoint(g->GetN(), x, F.GetParameter(NL + 1) * hps->Interpolate(x));
      }
      g->SetLineColor(kBlue + 1); g->SetLineWidth(3); g->SetLineStyle(2); g->Draw("L same");
      lg->AddEntry(g, "1n continuum (S_{n}=8.176)", "l");
   }
   lg->AddEntry(gb, "flat background", "l");
   auto *ls = new TLine(8.176, 0, 8.176, h->GetMaximum() * 0.55);
   ls->SetLineColor(kMagenta + 1); ls->SetLineWidth(2); ls->Draw();
   lg->AddEntry(ls, "S_{n} = 8.176", "l");
   lg->Draw();
   c->SaveAs(here + "/plots/fit_ps_" + tag + (usePS ? "" : "_noPS") + ".png");
}
