/// @file psa_sweep_ridge_C14.C
/// @brief Does the proton KE excess move when the hit content entering the fit is changed?
///
/// The +7-19 % excess above KE 8.5 MeV has survived beam energy, drift velocity, field scale,
/// vertex z, the chi2 cut, the acceptance, the estimator, the reconstruction chain in simulation,
/// and the clustering algorithm. Space charge is excluded quantitatively (lambda ~ 6.5e-12 C/m at
/// 2000 pps -> at most 0.23 mm of radial distortion, against a 0.5 mm hit sigma). What is left is
/// the point cloud itself: the simulation carries ~1.8x fewer pads per track than the data, and
/// above 8.5 MeV the helix no longer closes inside the detector, so the fit is a partial arc where
/// curvature, angle and momentum are nearly degenerate and a few stray hits weigh heavily.
///
/// psa_sweep_C14.sh reconstructs the same events at PSA thresholds 20/40/80 and with the
/// AtDirDeDxCleaner off. This macro reads the resulting Ex caches and reports, per variant:
///   * the elastic ridge offset from the kinematic line, averaged over theta_lab 34-58 deg
///     (where the bias lives) and over 66-78 deg (the control region, KE < 8.5 MeV, no bias)
///   * how many tracks survived, so a bias reduction bought by throwing away the sample is visible
///
/// Reading the result: if the bias falls as the threshold rises, stray hits are pulling the fit;
/// if it does not move while the statistics do, the hits themselves are mispositioned and the PSA
/// z-assignment is the next suspect.
///
///   root -b -q 'psa_sweep_ridge_C14.C()'

#include <tuple>

static double ps_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

static double ps_ex(double m1, double m2, double m3, double m4, double Eb, double thl, double Ke)
{
   double Et1 = Eb + m1, Et3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double a = (std::cos(thl) * ps_om2(s, m1 * m1, m2 * m2) * ps_om2(u, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                 (2 * m2 * m2) +
              s + u - m2 * m2;
   return a > 0 ? std::sqrt(a) - m4 : NAN;
}

static double ps_keFor(double m1, double m2, double Eb, double thl)
{
   double lo = 0.02, hi = 0.999 * Eb;
   auto g = [&](double k) { return ps_ex(m1, m2, m2, m1, Eb, thl, k); };
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

void psa_sweep_ridge_C14(Double_t Eb = 161.0)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   const double u = 931.49401, m1 = 14.003242 * u, m2 = 1.007825 * u;

   const int NV = 4;
   const char *cache[NV] = {"plots/proton_kin_pracmp_hdbscan.root", "plots/proton_kin_sweep_thr40.root",
                            "plots/proton_kin_sweep_thr80.root", "plots/proton_kin_sweep_noclean.root"};
   const char *lbl[NV] = {"thr 20 + clean", "thr 40 + clean", "thr 80 + clean", "thr 20, NO clean"};
   const int col[NV] = {kBlack, kAzure + 2, kRed + 1, kOrange + 7};
   const int mk[NV] = {20, 21, 22, 23};

   TTree *t[NV] = {nullptr, nullptr, nullptr, nullptr};
   double ntrk[NV] = {0, 0, 0, 0};
   for (int i = 0; i < NV; ++i) {
      TFile *f = TFile::Open(here + "/" + cache[i]);
      if (!f || f->IsZombie()) {
         printf("\033[1;33mmissing %s -- skipping this variant\033[0m\n", cache[i]);
         continue;
      }
      t[i] = (TTree *)f->Get("pk");
      if (t[i])
         ntrk[i] = t[i]->GetEntries();
   }

   auto *g = new TGraph[NV];
   int np[NV] = {0, 0, 0, 0};
   printf("\n===== elastic ridge offset from the kinematic line, Ebeam = %.0f =====\n", Eb);
   printf("  theta_lab | KE_line |");
   for (int i = 0; i < NV; ++i)
      printf(" %-16s |", lbl[i]);
   printf("\n");
   double sumBias[NV] = {0, 0, 0, 0}, sumCtl[NV] = {0, 0, 0, 0};
   int nBias[NV] = {0, 0, 0, 0}, nCtl[NV] = {0, 0, 0, 0};
   for (double th = 22; th < 82; th += 4) {
      double keL = ps_keFor(m1, m2, Eb, (th + 2.0) * TMath::DegToRad());
      if (!std::isfinite(keL))
         continue;
      printf("  %4.0f-%-4.0f | %7.2f |", th, th + 4, keL);
      for (int i = 0; i < NV; ++i) {
         if (!t[i]) {
            printf(" %16s |", "-");
            continue;
         }
         auto *h = new TH1D(TString::Format("hs%d_%d", i, (int)th), "", 140, std::max(0.0, 0.55 * keL), 1.55 * keL);
         t[i]->Draw(TString::Format("ke>>hs%d_%d", i, (int)th), TString::Format("theta>=%g&&theta<%g", th, th + 4),
                    "goff");
         h->SetDirectory(nullptr);
         if (h->Integral() < 40) {
            printf(" %16s |", "too few");
            delete h;
            continue;
         }
         h->Smooth(1);
         double pk = h->GetBinCenter(h->GetMaximumBin());
         double rel = 100.0 * (pk - keL) / keL;
         printf("  %7.2f %+6.1f%% |", pk, rel);
         g[i].SetPoint(np[i]++, th + 2, rel);
         double c = th + 2;
         if (c >= 34 && c <= 58) {
            sumBias[i] += rel;
            ++nBias[i];
         } else if (c >= 66 && c <= 78) {
            sumCtl[i] += rel;
            ++nCtl[i];
         }
         delete h;
      }
      printf("\n");
   }

   printf("\n  variant           | tracks  | mean bias 34-58 deg (KE>8.5) | mean 66-78 deg (control)\n");
   for (int i = 0; i < NV; ++i) {
      if (!t[i])
         continue;
      // Print the BIN COUNT beside each mean and refuse to average zero bins. Averaging an empty
      // set silently prints 0.0 %, which reads as "no bias" when it means "no data" -- exactly the
      // failure this table exists to expose.
      TString b = nBias[i] ? TString::Format("%+8.1f %% (%d bins)", sumBias[i] / nBias[i], nBias[i])
                           : TString("     n/a  (0 bins)");
      TString c = nCtl[i] ? TString::Format("%+8.1f %% (%d bins)", sumCtl[i] / nCtl[i], nCtl[i])
                          : TString("     n/a  (0 bins)");
      printf("  %-17s | %7.0f | %-27s | %s\n", lbl[i], ntrk[i], b.Data(), c.Data());
   }
   printf("\n  (a bias that falls only because the sample shrank is not a fix -- compare all three columns)\n");

   TCanvas *c1 = new TCanvas("c1", "psa sweep", 900, 620);
   bool first = true;
   auto *lg = new TLegend(0.50, 0.70, 0.89, 0.89);
   for (int i = 0; i < NV; ++i) {
      if (!t[i] || np[i] == 0)
         continue;
      g[i].SetMarkerStyle(mk[i]);
      g[i].SetMarkerColor(col[i]);
      g[i].SetLineColor(col[i]);
      g[i].SetLineWidth(2);
      g[i].SetMarkerSize(1.2);
      if (first) {
         g[i].SetTitle("elastic KE bias vs #theta_{lab}, same events;#theta_{lab} [deg];(ridge-line)/line [%]");
         g[i].GetYaxis()->SetRangeUser(-25, 30);
         g[i].GetXaxis()->SetLimits(20, 84);
         g[i].Draw("ALP");
         first = false;
      } else
         g[i].Draw("LP same");
      lg->AddEntry(&g[i], lbl[i], "lp");
   }
   auto *z = new TLine(20, 0, 84, 0);
   z->SetLineStyle(2);
   z->SetLineColor(kGray + 2);
   z->Draw();
   // the partial-arc threshold: KE 8.5 MeV is theta_lab ~ 62 deg
   auto *thr = new TLine(62, -25, 62, 30);
   thr->SetLineStyle(3);
   thr->SetLineColor(kGray + 2);
   thr->Draw();
   lg->Draw();
   TString png = here + "/plots/psa_sweep_ridge_C14.png";
   c1->SaveAs(png);
   printf("\nwrote %s\n\n", png.Data());
}
