/// @file plot_kt_2d.C
/// @brief KE vs theta for the (d,t) tritons, binned coarsely enough to actually read.
///
/// The first version used 90 x 80 bins over 70 deg x 12 MeV -- 0.78 deg x 0.15 MeV, about 4
/// tracks per bin for matFX-off and 2 for CATIMA, which renders as speckle rather than a ridge.
/// Default here is 44 x 34 (1.6 deg x 0.35 MeV), roughly 5x the occupancy per bin.
///
/// Loci are exact relativistic two-body solutions in the LAB frame (total energy E_beam + m_d,
/// target at rest -- NOT sqrt(s)). Two triton energies are allowed below the ~62 deg cap; the
/// analysis lives on the low one, so both branches are drawn, solid low / dashed high.
///
///   root -b -q 'plot_kt_2d.C(44, 34)'          // coarser: (30, 24)   finer: (60, 48)

#include <cmath>
#include <vector>

namespace {
const double u = 931.49410242;
const double m_C16 = 16.0147013 * u, m_d = 2.0135532 * u;
const double m_t = 3.01550072 * u, m_C15 = 15.0105993 * u;
const double Ebeam = 184.17;

void roots(double thetaRad, double Ex, double &kLo, double &kHi)
{
   kLo = kHi = -1;
   const double m4 = m_C15 + Ex;
   const double E1 = m_C16 + Ebeam, p1 = std::sqrt(E1 * E1 - m_C16 * m_C16);
   const double Etot = E1 + m_d, c = std::cos(thetaRad);
   auto f = [&](double p3) {
      return std::sqrt(p3 * p3 + m_t * m_t) +
             std::sqrt(p1 * p1 + p3 * p3 - 2 * p1 * p3 * c + m4 * m4) - Etot;
   };
   std::vector<double> r;
   double pPrev = 0, fPrev = f(0);
   const int NS = 20000;
   const double pMax = 1.2 * Etot;
   for (int i = 1; i <= NS; ++i) {
      const double p = pMax * i / NS, fv = f(p);
      if ((fPrev < 0) != (fv < 0)) {
         double lo = pPrev, hi = p;
         for (int it = 0; it < 80; ++it) {
            const double mid = 0.5 * (lo + hi);
            if ((f(lo) < 0) == (f(mid) < 0)) lo = mid; else hi = mid;
         }
         r.push_back(0.5 * (lo + hi));
      }
      pPrev = p; fPrev = fv;
   }
   if (r.empty()) return;
   std::sort(r.begin(), r.end());
   auto ke = [&](double p) { return std::sqrt(p * p + m_t * m_t) - m_t; };
   kLo = ke(r.front()); kHi = ke(r.back());
}

void drawLoci(double thMax, double keMax, bool labels)
{
   const double lv[] = {0.0, 0.740, 3.103, 4.220};
   for (int i = 0; i < 4; ++i) {
      auto *gLo = new TGraph(), *gHi = new TGraph();
      for (double th = 0.5; th <= 70; th += 0.5) {
         double a, b;
         roots(th * TMath::DegToRad(), lv[i], a, b);
         if (a > 0 && th <= thMax && a <= keMax) gLo->SetPoint(gLo->GetN(), th, a);
         if (b > 0 && th <= thMax && b <= keMax) gHi->SetPoint(gHi->GetN(), th, b);
      }
      gLo->SetLineColor(kBlack); gLo->SetLineWidth(3);
      gHi->SetLineColor(kGray + 3); gHi->SetLineWidth(2); gHi->SetLineStyle(2);
      if (gLo->GetN() > 1) gLo->Draw("L same");
      if (gHi->GetN() > 1) gHi->Draw("L same");
      if (labels && gLo->GetN() > 1) {
         double x, y; gLo->GetPoint(gLo->GetN() - 1, x, y);
         auto *t = new TLatex(x + 0.8, y + 0.1, Form("%.2f", lv[i]));
         t->SetTextSize(0.034); t->SetTextColor(kBlack); t->Draw();
      }
   }
}
} // namespace

void plot_kt_2d(int nTh = 44, int nKe = 34, double thMax = 70, double keMax = 12,
                TString out = "plots/kt_2d_rebinned.png")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(60);
   gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);

   auto *fN = TFile::Open("/mnt/f/a1975/caches/dt_kin_catima.root");
   auto *fO = TFile::Open("/mnt/f/a1975/caches/dt_kin_dv1104.root");
   auto *tN = (TTree *)fN->Get("pk");
   auto *tO = (TTree *)fO->Get("pk");

   struct Arm { TTree *t; const char *cut; const char *title; };
   std::vector<Arm> arms = {
      {tO, "chi2ndf<1e9",  "matFX OFF  (current production)"},
      {tN, "chi2ndf<1e9",  "CATIMA, matFX ON  (converged)"},
      {tN, "chi2ndf>=1e9", "CATIMA  (the collapsed fits)"}};

   auto *c = new TCanvas("c", "", 1900, 620);
   c->Divide(3, 1);
   printf("\nbinning: %d x %d over (0-%g deg, 0-%g MeV) = %.2f deg x %.2f MeV per bin\n",
          nTh, nKe, thMax, keMax, thMax / nTh, keMax / nKe);

   for (size_t i = 0; i < arms.size(); ++i) {
      c->cd(i + 1);
      gPad->SetLogz(); gPad->SetRightMargin(0.14);
      gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13);
      auto *h = new TH2D(Form("h%zu", i), Form("%s;#theta_{lab}  (deg);triton KE  (MeV)", arms[i].title),
                         nTh, 0, thMax, nKe, 0, keMax);
      arms[i].t->Draw(Form("ke:theta>>h%zu", i), arms[i].cut, "goff");
      h->GetXaxis()->SetTitleSize(0.050); h->GetYaxis()->SetTitleSize(0.050);
      h->GetXaxis()->SetLabelSize(0.045); h->GetYaxis()->SetLabelSize(0.045);
      h->GetXaxis()->SetTitleOffset(1.10); h->GetYaxis()->SetTitleOffset(1.15);
      h->Draw("colz");
      drawLoci(thMax, keMax, i == 0);
      printf("  %-34s entries in view: %8.0f   max/bin: %5.0f\n",
             arms[i].title, h->Integral(), h->GetMaximum());
   }
   c->SaveAs(out);
   printf("wrote %s\n\n", out.Data());
}
