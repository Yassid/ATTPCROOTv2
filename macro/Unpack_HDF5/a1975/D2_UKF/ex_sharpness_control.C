/// @file ex_sharpness_control.C
/// @brief Is the sharper CATIMA Ex spectrum the FITTER, or just a harsher selection?
///
/// The CATIMA production keeps 15361 of its 37145 tracks; the rest collapsed (ndf <= 0) and are
/// cut. The matFX-off production keeps 35904 of 36029. So the two spectra are built from
/// populations of very different purity, and a narrower peak could be selection rather than
/// physics -- exactly the confound that made an earlier all-pairs comparison claim a 42.5%
/// resolution win that was entirely a fit-failure tail.
///
/// CONTROL: cut the matFX-off sample down to the SAME SIZE on fit quality (tightest chi2/ndf),
/// and on track length (highest ncl), and refit the same peak. If those sharpen as much as
/// CATIMA does, the fitter is not the cause.
///
/// The peak fitted is the strong one near 3 MeV (15C 3.103). Gaussian + linear background,
/// Poisson likelihood, identical window and identical binning for every arm.
///
///   root -b -q ex_sharpness_control.C

#include <cstdio>

namespace {

struct Res { double n, mu, sig, dmu, dsig; };

Res fitPeak(TTree *t, const char *cut, const char *name, TVirtualPad *pad, int col, const char *title)
{
   auto *h = new TH1D(name, Form("%s;E_{x}(^{15}C)  (MeV);counts / 100 keV", title), 110, -3, 8);
   t->Draw(Form("ex>>%s", name), cut, "goff");
   h->SetLineColor(col);
   h->SetLineWidth(2);

   auto *f = new TF1(Form("f_%s", name), "gaus(0)+pol1(3)", 1.6, 4.4);
   const int b = h->FindBin(2.95);
   f->SetParameters(h->GetBinContent(b), 2.95, 0.30, h->GetBinContent(h->FindBin(4.3)), 0.0);
   f->SetParLimits(1, 2.3, 3.6);    // position: bracket the 3.103 level
   f->SetParLimits(2, 0.06, 1.00);  // width kept physical so it cannot eat the background
   h->Fit(f, "QRL");

   pad->cd();
   h->GetXaxis()->SetRangeUser(-1.5, 6.0);
   h->Draw("hist");
   f->Draw("same");
   auto *tx = new TLatex(0.14, 0.84, Form("N = %.0f", (double)h->GetEntries()));
   tx->SetNDC(); tx->SetTextSize(0.045); tx->Draw();
   auto *tx2 = new TLatex(0.14, 0.78, Form("#mu = %.3f #pm %.3f", f->GetParameter(1), f->GetParError(1)));
   tx2->SetNDC(); tx2->SetTextSize(0.045); tx2->Draw();
   auto *tx3 = new TLatex(0.14, 0.72, Form("#sigma = %.3f #pm %.3f", fabs(f->GetParameter(2)), f->GetParError(2)));
   tx3->SetNDC(); tx3->SetTextSize(0.045); tx3->Draw();

   return {(double)h->GetEntries(), f->GetParameter(1), fabs(f->GetParameter(2)),
           f->GetParError(1), f->GetParError(2)};
}

// threshold on `var` that leaves approximately `want` entries passing `base`
double thresholdFor(TTree *t, const char *base, const char *var, long want, bool keepSmall)
{
   double lo = 0, hi = 1e6;
   if (!keepSmall) { lo = 0; hi = 1e6; }
   for (int i = 0; i < 60; ++i) {
      const double mid = 0.5 * (lo + hi);
      const long n = t->GetEntries(Form("%s && %s %s %g", base, var, keepSmall ? "<" : ">", mid));
      if (keepSmall) { if (n > want) hi = mid; else lo = mid; }
      else           { if (n > want) lo = mid; else hi = mid; }
   }
   return 0.5 * (lo + hi);
}

} // namespace

void ex_sharpness_control()
{
   gStyle->SetOptStat(0);
   auto *fN = TFile::Open("/mnt/f/a1975/caches/dt_kin_catima.root");
   auto *fO = TFile::Open("/mnt/f/a1975/caches/dt_kin_dv1104.root");
   auto *tN = (TTree *)fN->Get("pk");
   auto *tO = (TTree *)fO->Get("pk");
   const char *conv = "chi2ndf<1e9";

   const long nCat = tN->GetEntries(conv);
   printf("\nCATIMA converged      : %ld\n", nCat);
   printf("matFX-off converged   : %lld\n", tO->GetEntries(conv));

   const double chiCut = thresholdFor(tO, conv, "chi2ndf", nCat, true);
   const double nclCut = thresholdFor(tO, conv, "ncl", nCat, false);
   printf("matched cuts on the matFX-off sample:  chi2ndf < %.4g   |   ncl > %.0f\n\n", chiCut, nclCut);

   auto *c = new TCanvas("c", "", 1600, 1200);
   c->Divide(2, 2);

   Res a = fitPeak(tO, conv, "hAll", c->cd(1), kAzure + 2, "matFX off - ALL converged");
   Res b = fitPeak(tN, conv, "hCat", c->cd(2), kOrange + 8, "CATIMA - all converged");
   Res d = fitPeak(tO, Form("%s && chi2ndf<%g", conv, chiCut), "hChi", c->cd(3), kAzure + 2,
                   "matFX off - size-matched on chi2/ndf");
   Res e = fitPeak(tO, Form("%s && ncl>%g", conv, nclCut), "hNcl", c->cd(4), kAzure + 2,
                   "matFX off - size-matched on ncl");

   c->SaveAs("plots/ex_sharpness_control.png");

   printf("=== 3 MeV peak, identical window and binning ===\n");
   printf("%-38s %8s %14s %14s\n", "sample", "N", "mu (MeV)", "sigma (MeV)");
   auto row = [](const char *s, Res r) {
      printf("%-38s %8.0f  %6.3f+-%.3f  %6.3f+-%.3f\n", s, r.n, r.mu, r.dmu, r.sig, r.dsig);
   };
   row("matFX off, all converged", a);
   row("CATIMA, all converged", b);
   row("matFX off, size-matched on chi2/ndf", d);
   row("matFX off, size-matched on ncl", e);
   printf("\nIf the size-matched matFX-off rows reach CATIMA's sigma, the sharpening is\n"
          "SELECTION. If they stay wide, it is the fitter.\n\n");
}
