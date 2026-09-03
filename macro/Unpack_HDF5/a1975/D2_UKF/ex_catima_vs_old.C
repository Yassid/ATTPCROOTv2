/// @file ex_catima_vs_old.C
/// @brief 15C excitation-energy spectrum: CATIMA material-effects production vs the old one.
///
/// The two caches are built by scripts identical in every selection, so they differ only by the
/// fitter configuration and can be overlaid directly.
///
/// ONE CUT MATTERS AND IT IS NOT OPTIONAL. ex_dt_a1975.C writes chi2ndf = 1e9 when ndf <= 0, and
/// its own chi2 cut is `c2n > chi2Cut` with chi2Cut = 1e9, which 1e9 does NOT satisfy -- so
/// COLLAPSED FITS ARE IN THE CACHE. They carry kinematics like any other track. Requiring
/// chi2ndf < 1e9 is the convergence cut in cache space and has to be applied to both arms.
///
/// Peaks are the 15C ground state (1/2+) and the 0.740 MeV first excited state (5/2+), so the
/// separation is a KNOWN 0.740 MeV and is the figure of merit: every chain so far has undershot
/// it by 0.3-0.7 MeV.
///
///   root -b -q ex_catima_vs_old.C

#include <cstdio>

namespace {

struct Fit {
   double p0, s0, p1, s1, sep, n;
};

Fit fitPair(TTree *t, const char *cut, const char *name, int col, TVirtualPad *pad)
{
   TH1D *h = new TH1D(name, ";E_{x}(^{15}C) [MeV];counts / 100 keV", 60, -1.5, 4.5);
   t->Draw(Form("ex>>%s", name), cut, "goff");
   h->SetLineColor(col);
   h->SetLineWidth(2);

   // two gaussians on a linear background; start the pair at the known 0 / 0.740
   TF1 *f = new TF1(Form("f_%s", name), "gaus(0)+gaus(3)+pol1(6)", -1.0, 2.5);
   const double peak = h->GetMaximum();
   f->SetParameters(peak, 0.0, 0.35, 0.5 * peak, 0.74, 0.35, 0.1 * peak, 0.0);
   f->SetParLimits(1, -0.8, 0.8);  // g.s. position
   f->SetParLimits(2, 0.10, 1.20); // widths kept physical so the fit cannot swallow the pair
   f->SetParLimits(4, 0.10, 2.20); // excited state
   f->SetParLimits(5, 0.10, 1.20);
   h->Fit(f, "QRL"); // Poisson likelihood: the tails here are thin

   pad->cd();
   h->Draw("hist");
   f->Draw("same");

   Fit r;
   r.p0 = f->GetParameter(1);
   r.s0 = fabs(f->GetParameter(2));
   r.p1 = f->GetParameter(4);
   r.s1 = fabs(f->GetParameter(5));
   r.sep = r.p1 - r.p0;
   r.n = h->GetEntries();
   return r;
}

void report(const char *tag, const Fit &f)
{
   printf("  %-26s N=%-7.0f  gs=%+6.3f (sig %5.3f)   ex1=%+6.3f (sig %5.3f)   SEP=%6.3f MeV  [true 0.740]\n", tag, f.n,
          f.p0, f.s0, f.p1, f.s1, f.sep);
}

} // namespace

void ex_catima_vs_old()
{
   auto *fN = TFile::Open("/mnt/f/a1975/caches/dt_kin_catima.root");
   auto *fO = TFile::Open("/mnt/f/a1975/caches/dt_kin_dv1104.root");
   auto *tN = (TTree *)fN->Get("pk");
   auto *tO = (TTree *)fO->Get("pk");

   printf("\nraw entries: CATIMA %lld,  old %lld\n", tN->GetEntries(), tO->GetEntries());
   printf("converged  : CATIMA %lld,  old %lld   (chi2ndf < 1e9)\n", tN->GetEntries("chi2ndf<1e9"),
          tO->GetEntries("chi2ndf<1e9"));

   auto *c = new TCanvas("c", "", 1400, 900);
   c->Divide(2, 2);

   printf("\n=== convergence cut only ===\n");
   report("CATIMA  matFX on", fitPair(tN, "chi2ndf<1e9", "hN0", kRed + 1, c->cd(1)));
   report("old     matFX off", fitPair(tO, "chi2ndf<1e9", "hO0", kBlue + 1, c->cd(2)));

   const char *q = "chi2ndf<1e9 && chi2ndf<10 && ncl>30";
   printf("\n=== plus a quality cut (chi2/ndf < 10, ncl > 30) ===\n");
   report("CATIMA  matFX on", fitPair(tN, q, "hN1", kRed + 1, c->cd(3)));
   report("old     matFX off", fitPair(tO, q, "hO1", kBlue + 1, c->cd(4)));

   c->SaveAs("plots/ex_catima_vs_old.png");
   printf("\nwrote plots/ex_catima_vs_old.png\n\n");
}
