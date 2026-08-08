/// @file pra_ridge_C14.C
/// @brief triplclust vs HDBSCAN: elastic ridge against the kinematic line, same events.
///
/// The paper (Ayyad et al., EPJ A 59:294) used triplclust; the current a1954 production uses
/// HDBSCAN. Clustering sets the points the curvature is fitted through, so it is the leading
/// remaining candidate for the +7-19 % proton KE excess seen for theta_lab < 62 deg -- the one
/// thing not yet excluded after beam energy, drift velocity, field scale, vertex z and the
/// simulation all came back clean.
///
/// Reads the two caches written by pra_cmp_finish_C14.sh (same run, same events, same PSA,
/// cleaner, par file and UKF settings -- only praType differs) and prints the ridge offset that
/// kine_plane_C14.C measured on the production data.
///
///   root -b -q 'pra_ridge_C14.C()'

#include <tuple>

static double pr_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

static double pr_ex(double m1, double m2, double m3, double m4, double Eb, double thl, double Ke)
{
   double Et1 = Eb + m1, Et3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double a = (std::cos(thl) * pr_om2(s, m1 * m1, m2 * m2) * pr_om2(u, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                 (2 * m2 * m2) +
              s + u - m2 * m2;
   return a > 0 ? std::sqrt(a) - m4 : NAN;
}

static double pr_keFor(double m1, double m2, double Eb, double thl, double exW)
{
   double lo = 0.02, hi = 0.999 * Eb;
   auto g = [&](double k) { return pr_ex(m1, m2, m2, m1, Eb, thl, k) - exW; };
   double flo = g(lo), fhi = g(hi);
   if (!std::isfinite(flo) || !std::isfinite(fhi) || flo * fhi > 0)
      return NAN;
   for (int i = 0; i < 120; ++i) {
      double mid = 0.5 * (lo + hi), fm = g(mid);
      if (!std::isfinite(fm))
         return NAN;
      if (flo * fm <= 0)
         hi = mid;
      else {
         lo = mid;
         flo = fm;
      }
   }
   return 0.5 * (lo + hi);
}

void pra_ridge_C14(Double_t Eb = 161.0)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   const double u = 931.49401, m1 = 14.003242 * u, m2 = 1.007825 * u;

   const int NS = 2;
   const char *cache[NS] = {"plots/proton_kin_pracmp_tc.root", "plots/proton_kin_pracmp_hdbscan.root"};
   const char *lbl[NS] = {"triplclust", "HDBSCAN"};
   const int col[NS] = {kGreen + 3, kViolet + 1};

   TTree *t[NS];
   TFile *f[NS];
   for (int i = 0; i < NS; ++i) {
      f[i] = TFile::Open(here + "/" + cache[i]);
      if (!f[i] || f[i]->IsZombie()) {
         printf("\033[1;31mmissing %s -- run pra_cmp_finish_C14.sh first\033[0m\n", cache[i]);
         return;
      }
      t[i] = (TTree *)f[i]->Get("pk");
      printf("%-12s : %lld tracks\n", lbl[i], t[i]->GetEntries());
   }

   printf("\n===== elastic ridge vs the kinematic line (Ebeam = %.0f) =====\n", Eb);
   printf("  theta_lab | KE_line |  triplclust  diff   rel   |  HDBSCAN   diff   rel\n");
   auto *g0 = new TGraph(), *g1 = new TGraph();
   int n0 = 0, n1 = 0;
   for (double th = 22; th < 84; th += 4) {
      double keL = pr_keFor(m1, m2, Eb, (th + 2.0) * TMath::DegToRad(), 0.0);
      if (!std::isfinite(keL))
         continue;
      printf("  %4.0f-%-4.0f | %7.2f |", th, th + 4, keL);
      for (int i = 0; i < NS; ++i) {
         auto *h = new TH1D(TString::Format("hp%d_%d", i, (int)th), "", 140, std::max(0.0, 0.55 * keL), 1.55 * keL);
         t[i]->Draw(TString::Format("ke>>hp%d_%d", i, (int)th), TString::Format("theta>=%g&&theta<%g", th, th + 4),
                    "goff");
         h->SetDirectory(nullptr);
         if (h->Integral() < 40) {
            printf(" %26s |", "too few");
            delete h;
            continue;
         }
         h->Smooth(1);
         double pk = h->GetBinCenter(h->GetMaximumBin());
         printf(" %10.2f %+6.2f %+6.1f%% |", pk, pk - keL, 100 * (pk - keL) / keL);
         if (i == 0) {
            g0->SetPoint(n0++, th + 2, 100 * (pk - keL) / keL);
         } else {
            g1->SetPoint(n1++, th + 2, 100 * (pk - keL) / keL);
         }
         delete h;
      }
      printf("\n");
   }

   TCanvas *c1 = new TCanvas("c1", "pra ridge", 900, 620);
   g0->SetMarkerStyle(20);
   g0->SetMarkerColor(col[0]);
   g0->SetLineColor(col[0]);
   g0->SetLineWidth(2);
   g1->SetMarkerStyle(21);
   g1->SetMarkerColor(col[1]);
   g1->SetLineColor(col[1]);
   g1->SetLineWidth(2);
   g0->SetTitle("elastic KE bias vs #theta_{lab}, same events;#theta_{lab} [deg];(ridge - line)/line [%]");
   g0->GetYaxis()->SetRangeUser(-25, 30);
   g0->Draw("ALP");
   g1->Draw("LP same");
   auto *z = new TLine(20, 0, 86, 0);
   z->SetLineStyle(2);
   z->SetLineColor(kGray + 2);
   z->Draw();
   auto *lg = new TLegend(0.58, 0.74, 0.89, 0.89);
   lg->AddEntry(g0, lbl[0], "lp");
   lg->AddEntry(g1, lbl[1], "lp");
   lg->Draw();
   TString png = here + "/plots/pra_ridge_C14.png";
   c1->SaveAs(png);
   printf("\nwrote %s\n\n", png.Data());
}
