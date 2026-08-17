/// @file plot_ex_counts.C
/// @brief 15C excitation energy in COUNTS, both productions overlaid.
///
/// Counts rather than a normalised fraction is the honest comparison now: after the TGeo
/// navigator fix the two samples are essentially the same size (35601 CATIMA vs 35904 matFX
/// off), so no normalisation is needed and none is applied. Before the fix CATIMA had 15361 and
/// the plots HAD to be area-normalised, which hid the yield difference.
///
/// Poisson error bars are drawn so peak significance can be judged by eye.
///
///   root -b -q 'plot_ex_counts.C'                       // -3 to 8 MeV
///   root -b -q 'plot_ex_counts.C(-1,6,0.10)'            // zoom, 100 keV bins
///   root -b -q 'plot_ex_counts.C(-1,6,0.05)'            // 50 keV bins

#include <cstdio>

void plot_ex_counts(double lo = -3, double hi = 8, double binMeV = 0.10,
                    TString out = "plots/ex_counts.png")
{
   gStyle->SetOptStat(0);
   gStyle->SetPadTickX(1);
   gStyle->SetPadTickY(1);

   auto *fN = TFile::Open("/mnt/f/a1975/caches/dt_kin_catima.root");
   auto *fO = TFile::Open("/mnt/f/a1975/caches/dt_kin_dv1104.root");
   if (!fN || !fO || fN->IsZombie() || fO->IsZombie()) { printf("cannot open caches\n"); return; }
   auto *tN = (TTree *)fN->Get("pk");
   auto *tO = (TTree *)fO->Get("pk");
   const char *conv = "chi2ndf<1e9";

   const int nb = (int)std::lround((hi - lo) / binMeV);
   auto *hO = new TH1D("hO", Form(";E_{x}(^{15}C)  (MeV);counts / %.0f keV", binMeV * 1000), nb, lo, hi);
   auto *hN = new TH1D("hN", "", nb, lo, hi);
   tO->Draw("ex>>hO", conv, "goff");
   tN->Draw("ex>>hN", conv, "goff");

   hO->SetLineColor(kAzure + 2); hO->SetLineWidth(2); hO->SetMarkerColor(kAzure + 2);
   hN->SetLineColor(kOrange + 8); hN->SetLineWidth(2); hN->SetMarkerColor(kOrange + 8);
   hO->Sumw2(); hN->Sumw2();

   auto *c = new TCanvas("c", "", 1500, 850);
   c->SetLeftMargin(0.11); c->SetBottomMargin(0.12);
   hO->GetXaxis()->SetTitleSize(0.045); hO->GetYaxis()->SetTitleSize(0.045);
   hO->GetXaxis()->SetLabelSize(0.040); hO->GetYaxis()->SetLabelSize(0.040);
   hO->GetYaxis()->SetTitleOffset(1.15);
   hO->SetMaximum(1.28 * std::max(hO->GetMaximum(), hN->GetMaximum()));
   hO->SetMinimum(0);
   hO->Draw("E0");
   hN->Draw("E0 same");

   // known 15C levels
   const double lv[] = {0.0, 0.740, 3.103, 4.220, 4.657};
   for (double x : lv) {
      if (x < lo || x > hi) continue;
      auto *l = new TLine(x, 0, x, hO->GetMaximum());
      l->SetLineStyle(3); l->SetLineColor(kGray + 2); l->Draw();
      auto *t = new TLatex(x + 0.03, 1.20 * std::max(hO->GetMaximum(), hN->GetMaximum()) / 1.28,
                           Form("%.3f", x));
      t->SetTextSize(0.026); t->SetTextColor(kGray + 3); t->SetTextAngle(90); t->Draw();
   }

   auto *leg = new TLegend(0.60, 0.75, 0.89, 0.89);
   leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
   leg->AddEntry(hO, Form("matFX off  -  %.0f tracks", hO->GetEntries()), "le");
   leg->AddEntry(hN, Form("CATIMA     -  %.0f tracks", hN->GetEntries()), "le");
   leg->Draw();

   c->SaveAs(out);
   printf("\nwrote %s   (%d bins of %.0f keV)\n", out.Data(), nb, binMeV * 1000);
   printf("  matFX off : %8.0f entries, peak bin %5.0f counts\n", hO->GetEntries(), hO->GetMaximum());
   printf("  CATIMA    : %8.0f entries, peak bin %5.0f counts\n", hN->GetEntries(), hN->GetMaximum());

   // counts in a window around each known level, both arms, so the gain is a number not an impression
   printf("\n%-10s %14s %14s %10s\n", "level", "matFX off", "CATIMA", "ratio");
   for (double x : lv) {
      if (x < lo || x > hi) continue;
      const int b1 = hO->FindBin(x - 0.25), b2 = hO->FindBin(x + 0.25);
      const double a = hO->Integral(b1, b2), b = hN->Integral(b1, b2);
      printf("%-10.3f %14.0f %14.0f %10.2f\n", x, a, b, a > 0 ? b / a : 0.0);
   }
   printf("\n(window is +-250 keV; these are RAW counts, no background subtraction)\n\n");
}
