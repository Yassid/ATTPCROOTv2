/// @file globalfit_check_C14.C
/// @brief ROOT reference for the explorer's "global fit" panel: the same model, fitted by Minuit.
///
/// Levels (positions literature + shift, widths sig0 + dSig*(E - 6.094)), the 14N (p,n) blend, the
/// 1n phase space (hPS of phasespace_1n_C14.C) and a linear background; POISSON likelihood, only the
/// amplitudes float. Prints the binned spectrum as well, so the page's live Ex can be compared bin by
/// bin, and the likelihood written out explicitly rather than Minuit's MinFcnValue, whose offset
/// convention differs between versions.
///
///   root -b -q 'globalfit_check_C14.C("plots/proton_kin_cat5_s013.root")'
namespace gc {
int NL = 5;
double E[8] = {6.091, 6.728, 7.012, 7.341, 8.317}, SG[8];
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
} // namespace gc

void globalfit_check_C14(TString cache = "plots/proton_kin_cat5_s013.root", Double_t cmLo = 20,
                         Double_t cmHi = 140, Double_t exLo = 5.0, Double_t exHi = 10.4,
                         Double_t binW = 0.05, Double_t vzLo = 10, Double_t vzHi = 490,
                         Double_t chi2Cut = 5.0, Double_t sig0 = 0.132, Double_t dSig = 0.0123)
{
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   TFile *fd = TFile::Open(here + "/" + cache);
   TNtuple *t = fd && !fd->IsZombie() ? (TNtuple *)fd->Get("pk") : nullptr;
   if (!t) { printf("no cache %s\n", cache.Data()); return; }
   TFile *fp = TFile::Open(here + "/plots/phasespace_1n_C14.root");
   gc::hps = (TH1D *)fp->Get("hPS");          // NB "hPS": spec_states_C14.C asks for "hps" and gets nothing
   gc::hps->SetDirectory(nullptr);
   gc::hps->Scale(1.0 / gc::hps->Integral());
   for (int i = 0; i < gc::NL; ++i) gc::SG[i] = sig0 + dSig * (gc::E[i] - 6.094);

   const int nb = (int)std::lround((exHi - exLo) / binW);
   auto *h = new TH1D("hgc", "", nb, exLo, exHi);
   t->Draw("ex>>hgc", Form("chi2ndf<%g && thcm>=%g && thcm<=%g && vertexz>=%g && vertexz<=%g",
                          chi2Cut, cmLo, cmHi, vzLo, vzHi), "goff");
   printf("BINS");
   for (int b = 1; b <= nb; ++b) printf(" %.0f", h->GetBinContent(b));
   printf("\nN %.0f\n", h->Integral());

   const int NL = gc::NL, NPAR = NL + 4;
   auto *F = new TF1("Fgc", gc::model, exLo, exHi, NPAR);
   for (int i = 0; i < NL; ++i) {
      F->SetParameter(i, std::max(1.0, h->GetBinContent(h->FindBin(gc::E[i] + gc::shift)) * 0.5));
      F->SetParLimits(i, 0, 1e5);
   }
   F->SetParameter(NL, 20);                    F->SetParLimits(NL, 0, 1e5);
   F->SetParameter(NL + 1, h->Integral() * 0.05); F->SetParLimits(NL + 1, 0, 1e7);
   F->SetParameter(NL + 2, 2);                 F->SetParameter(NL + 3, 0);
   TFitResultPtr r = h->Fit(F, "RQNSL");
   for (int pass = 0; pass < 3 && ((int)r->Status() != 0 || r->CovMatrixStatus() < 2); ++pass) r = h->Fit(F, "RQNSL");
   printf("status %d cov %d\n", (int)r->Status(), r->CovMatrixStatus());

   double nll = 0, bc = 0, pear = 0; int nfree = 0;
   for (int b = 1; b <= nb; ++b) {
      double m = F->Eval(h->GetBinCenter(b)), n = h->GetBinContent(b);
      nll += m - (n > 0 ? n * std::log(m) : 0);
      bc += 2 * (m - n + (n > 0 ? n * std::log(n / m) : 0));
      if (m > 0) { pear += (n - m) * (n - m) / m; ++nfree; }
   }
   const double sq = std::sqrt(2 * TMath::Pi());
   for (int i = 0; i < NL; ++i)
      printf("LEVEL %.3f area %.2f err %.2f\n", gc::E[i], F->GetParameter(i) * gc::SG[i] * sq / binW,
             F->GetParError(i) * gc::SG[i] * sq / binW);
   printf("N14 area %.2f err %.2f\n", F->GetParameter(NL) * gc::sgN * sq / binW, F->GetParError(NL) * gc::sgN * sq / binW);
   double psInt = 0;
   for (int b = 1; b <= nb; ++b) psInt += F->GetParameter(NL + 1) * gc::hps->Interpolate(h->GetBinCenter(b));
   printf("PS counts %.2f\n", psInt);
   printf("BG %.4f + %.4f x  -> at lo %.3f at hi %.3f, counts %.2f\n", F->GetParameter(NL + 2), F->GetParameter(NL + 3),
          F->GetParameter(NL + 2) + F->GetParameter(NL + 3) * exLo, F->GetParameter(NL + 2) + F->GetParameter(NL + 3) * exHi,
          nb * (F->GetParameter(NL + 2) + F->GetParameter(NL + 3) * 0.5 * (exLo + exHi)));
   printf("NLL %.4f  BakerCousins %.4f  Pearson %.3f ndf %d chi2/ndf %.3f\n", nll, bc, pear, nfree - NPAR,
          pear / (nfree - NPAR));
}
