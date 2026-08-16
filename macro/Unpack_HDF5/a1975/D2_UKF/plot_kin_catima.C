/// @file plot_kin_catima.C
/// @brief (d,t) kinematics: CATIMA material-effects production vs the matFX-off production.
///
/// Writes a single PNG. The kinematic loci are EXACT relativistic two-body solutions in the LAB
/// frame -- total energy is E_beam + m_d with the target at rest, NOT sqrt(s); using the CM
/// energy there yields no roots at all. Below the ~62 deg cap TWO triton energies are allowed at
/// each angle, and the (d,t) analysis lives on the LOW one, so both branches are drawn.
///
///   root -b -q plot_kin_catima.C

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

void drawLoci(double thMax, double keMax)
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
      gLo->SetLineColor(kBlack); gLo->SetLineWidth(3); gLo->SetLineStyle(1);
      gHi->SetLineColor(kGray + 2); gHi->SetLineWidth(2); gHi->SetLineStyle(2);
      if (gLo->GetN() > 1) gLo->Draw("L same");
      if (gHi->GetN() > 1) gHi->Draw("L same");
      if (gLo->GetN() > 1) {
         double x, y; gLo->GetPoint(gLo->GetN() - 1, x, y);
         auto *t = new TLatex(x + 1, y, Form("%.2f", lv[i]));
         t->SetTextSize(0.030); t->SetTextColor(kBlack); t->Draw();
      }
   }
}

TH2D *kt(TTree *t, const char *name, const char *cut, double thMax, double keMax)
{
   auto *h = new TH2D(name, "", 90, 0, thMax, 80, 0, keMax);
   t->Draw(Form("ke:theta>>%s", name), cut, "goff");
   h->GetXaxis()->SetTitle("#theta_{lab}  (deg)");
   h->GetYaxis()->SetTitle("triton KE  (MeV)");
   h->GetXaxis()->SetTitleSize(0.045); h->GetYaxis()->SetTitleSize(0.045);
   h->GetXaxis()->SetLabelSize(0.040); h->GetYaxis()->SetLabelSize(0.040);
   return h;
}
} // namespace

void plot_kin_catima(TString out = "plots/kinematics_catima.png")
{
   gStyle->SetOptStat(0);
   gStyle->SetPalette(kBird);
   gStyle->SetNumberContours(64);

   auto *fN = TFile::Open("/mnt/f/a1975/caches/dt_kin_catima.root");
   auto *fO = TFile::Open("/mnt/f/a1975/caches/dt_kin_dv1104.root");
   if (!fN || !fO || fN->IsZombie() || fO->IsZombie()) { printf("cannot open caches\n"); return; }
   auto *tN = (TTree *)fN->Get("pk");
   auto *tO = (TTree *)fO->Get("pk");
   const char *conv = "chi2ndf<1e9";

   auto *c = new TCanvas("c", "", 1600, 1200);
   c->Divide(2, 2);

   // ── 1/2: KE vs theta, zoomed to where the physics is ──
   c->cd(1); gPad->SetLogz(); gPad->SetRightMargin(0.13); gPad->SetLeftMargin(0.12);
   auto *hO = kt(tO, "hO", conv, 70, 12);
   hO->SetTitle("matFX OFF  (current production)");
   hO->Draw("colz"); drawLoci(70, 12);

   c->cd(2); gPad->SetLogz(); gPad->SetRightMargin(0.13); gPad->SetLeftMargin(0.12);
   auto *hN = kt(tN, "hN", conv, 70, 12);
   hN->SetTitle("CATIMA, matFX ON");
   hN->Draw("colz"); drawLoci(70, 12);

   // ── 3: Ex spectra, area-normalised ──
   c->cd(3); gPad->SetLeftMargin(0.12);
   auto *eO = new TH1D("eO", ";E_{x}(^{15}C)  (MeV);fraction / 100 keV", 110, -3, 8);
   auto *eN = new TH1D("eN", "", 110, -3, 8);
   tO->Draw("ex>>eO", conv, "goff");
   tN->Draw("ex>>eN", conv, "goff");
   if (eO->Integral() > 0) eO->Scale(1.0 / eO->Integral());
   if (eN->Integral() > 0) eN->Scale(1.0 / eN->Integral());
   eO->SetLineColor(kAzure + 2); eO->SetLineWidth(3);
   eN->SetLineColor(kOrange + 8); eN->SetLineWidth(3);
   eO->SetMaximum(1.25 * std::max(eO->GetMaximum(), eN->GetMaximum()));
   eO->GetXaxis()->SetTitleSize(0.045); eO->GetYaxis()->SetTitleSize(0.045);
   eO->SetTitle("^{15}C excitation energy");
   eO->Draw("hist"); eN->Draw("hist same");
   for (double x : {0.0, 0.740, 3.103, 4.220}) {
      auto *l = new TLine(x, 0, x, eO->GetMaximum());
      l->SetLineStyle(3); l->SetLineColor(kGray + 2); l->Draw();
   }
   auto *leg = new TLegend(0.56, 0.74, 0.88, 0.88);
   leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.034);
   leg->AddEntry(eO, Form("matFX off  (%.0f)", (double)tO->GetEntries(conv)), "l");
   leg->AddEntry(eN, Form("CATIMA     (%.0f)", (double)tN->GetEntries(conv)), "l");
   leg->Draw();

   // ── 4: per-run collapse vs vertex z ──
   c->cd(4); gPad->SetLeftMargin(0.12);
   int runs[] = {16, 17, 18, 19, 20, 21, 22, 23, 26, 27, 31, 32, 34, 36, 37, 38, 39, 40, 41, 42,
                 43, 44, 46, 48, 57, 58, 76, 77, 78, 79, 80, 82, 83, 84, 85, 86, 87, 88, 89, 91,
                 92, 95, 96, 97, 98, 102, 103};
   auto *gGood = new TGraph(), *gBad = new TGraph();
   std::vector<double> X, Y;
   for (int r : runs) {
      long n = tN->GetEntries(Form("run==%d", r));
      if (!n) continue;
      long b = tN->GetEntries(Form("run==%d && chi2ndf>=1e9", r));
      auto *hz = new TH1D("hz", "", 100, 0, 1000);
      tN->Draw("vertexz>>hz", Form("run==%d && chi2ndf<1e9", r), "goff");
      const double vz = hz->GetMean(), pct = 100.0 * b / n;
      ((pct < 20) ? gGood : gBad)->SetPoint(((pct < 20) ? gGood : gBad)->GetN(), vz, pct);
      X.push_back(vz); Y.push_back(pct);
      delete hz;
   }
   auto *frame = new TH2D("frame", "per-run fit collapse;mean converged vertex z  (mm);fits collapsed  (%)",
                          10, 300, 540, 10, -5, 105);
   frame->GetXaxis()->SetTitleSize(0.045); frame->GetYaxis()->SetTitleSize(0.045);
   frame->Draw();
   gBad->SetMarkerStyle(20); gBad->SetMarkerSize(1.4); gBad->SetMarkerColor(kRed + 1);
   gGood->SetMarkerStyle(20); gGood->SetMarkerSize(1.4); gGood->SetMarkerColor(kAzure + 2);
   gBad->Draw("P same"); gGood->Draw("P same");
   const int n = X.size();
   double sx = 0, sy = 0, sxx = 0, sxy = 0, syy = 0;
   for (int i = 0; i < n; ++i) { sx += X[i]; sy += Y[i]; sxx += X[i] * X[i];
      sxy += X[i] * Y[i]; syy += Y[i] * Y[i]; }
   const double sl = (n * sxy - sx * sy) / (n * sxx - sx * sx), ic = (sy - sl * sx) / n;
   const double rr = (n * sxy - sx * sy) / std::sqrt((n * sxx - sx * sx) * (n * syy - sy * sy));
   auto *tl = new TLine(300, ic + sl * 300, 540, ic + sl * 540);
   tl->SetLineStyle(2); tl->SetLineColor(kGray + 2); tl->SetLineWidth(2); tl->Draw();
   auto *lab = new TLatex(0.60, 0.20, Form("r = %+.2f,  n = %d", rr, n));
   lab->SetNDC(); lab->SetTextSize(0.036); lab->Draw();

   c->SaveAs(out);
   printf("\nwrote %s\n", out.Data());
   printf("  matFX off : %lld converged of %lld\n", tO->GetEntries(conv), tO->GetEntries());
   printf("  CATIMA    : %lld converged of %lld\n", tN->GetEntries(conv), tN->GetEntries());
   printf("  collapse correlation with vertex z: r = %+.3f (n = %d)\n", rr, n);
}
