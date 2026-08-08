/// @file dv_emulate_C14.C
/// @brief Emulate a drift-velocity change on already-fitted tracks, without re-reconstructing.
///
/// The drift velocity converts time buckets to z. If the reconstruction used dv_used but the true
/// value is dv_true = k * dv_used, then every reconstructed z is short by the factor k while the
/// TRANSVERSE geometry -- and therefore the curvature and p_T -- is untouched:
///
///     tan(theta_true) = r_T / (k * z_rec) = tan(theta_rec) / k
///     p_T             = p_rec * sin(theta_rec)        (invariant)
///     p_true          = p_T / sin(theta_true)
///
/// So a dv error maps (KE, theta) to a new (KE, theta) with no refit, to first order: it ignores
/// that the fit itself would have converged slightly differently, but it captures the dominant
/// geometric effect exactly. That is enough to ask which k flattens the elastic E_x locus, which
/// costs seconds instead of the ~2 h a full re-reco per dv point would take.
///
/// Flatness metric: the spread of the elastic E_x peak position across theta_lab bands spanning
/// 24-76 deg. The bands are chosen in theta_LAB because that is what the drift error acts on.
///
///   root -b -q 'dv_emulate_C14.C()'

#include <tuple>

static double dv_om2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

static double dv_ex(double m1, double m2, double m3, double m4, double Eb, double thl, double Ke)
{
   double Et1 = Eb + m1, Et3 = Ke + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double a = (std::cos(thl) * dv_om2(s, m1 * m1, m2 * m2) * dv_om2(u, m2 * m2, m3 * m3) -
               (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                 (2 * m2 * m2) +
              s + u - m2 * m2;
   return a > 0 ? std::sqrt(a) - m4 : NAN;
}

void dv_emulate_C14(Double_t dvUsed = 1.30, Double_t Eb = 161.0, Double_t kLo = 0.85, Double_t kHi = 1.35,
                    Double_t dk = 0.025)
{
   gStyle->SetOptStat(0);
   TString here = gSystem->DirName(gInterpreter->GetCurrentMacroName());
   const double u = 931.49401, m1 = 14.003242 * u, m2 = 1.007825 * u, m_p = m2;

   const int NS = 2;
   const char *file[NS] = {"plots/proton_kin_300_ukf_nc.root", "plots/proton_kin_300gfx_nc.root"};
   const char *lbl[NS] = {"UKF", "GENFIT"};
   const int col[NS] = {kAzure + 2, kRed + 1};

   const int NB = 9;
   int blo[NB] = {24, 30, 36, 42, 48, 54, 60, 66, 72};
   const int bw = 6;

   auto *gSpread = new TGraph[NS];
   auto *gMuAvg = new TGraph[NS];
   int npg[NS] = {0, 0};
   double bestK[NS] = {1, 1}, bestS[NS] = {1e9, 1e9};

   for (int i = 0; i < NS; ++i) {
      TFile *f = TFile::Open(here + "/" + file[i]);
      if (!f || f->IsZombie()) {
         printf("\033[1;31mmissing %s\033[0m\n", file[i]);
         return;
      }
      TTree *t = (TTree *)f->Get("pk");
      float ke, th;
      t->SetBranchAddress("ke", &ke);
      t->SetBranchAddress("theta", &th);
      Long64_t N = t->GetEntries();
      // cache the tracks once
      std::vector<float> vKE, vTh;
      vKE.reserve(N);
      vTh.reserve(N);
      for (Long64_t e = 0; e < N; ++e) {
         t->GetEntry(e);
         if (ke > 0.2 && ke < 200 && th > 5 && th < 100) {
            vKE.push_back(ke);
            vTh.push_back(th);
         }
      }
      f->Close();
      printf("\n===== %s : %zu tracks, dv_used = %.3f, Ebeam = %.1f =====\n", lbl[i], vKE.size(), dvUsed, Eb);
      printf("     k    dv_true |");
      for (int b = 0; b < NB; ++b)
         printf(" %2d-%2d", blo[b], blo[b] + bw);
      printf(" | spread\n");

      for (double k = kLo; k <= kHi + 1e-9; k += dk) {
         std::vector<TH1D *> h(NB);
         for (int b = 0; b < NB; ++b) {
            h[b] = new TH1D(TString::Format("hh%d_%d_%d", i, b, (int)std::lround(k * 1000)), "", 140, -4, 3);
            h[b]->SetDirectory(nullptr);
         }
         for (size_t e = 0; e < vKE.size(); ++e) {
            double thR = vTh[e] * TMath::DegToRad();
            // geometric remap
            double thT = std::atan2(std::sin(thR), k * std::cos(thR));
            double pR = std::sqrt(vKE[e] * (vKE[e] + 2 * m_p));
            double pT = pR * std::sin(thR);
            double sT = std::sin(thT);
            if (sT < 1e-6)
               continue;
            double pN = pT / sT;
            double keN = std::sqrt(pN * pN + m_p * m_p) - m_p;
            double thD = thT * TMath::RadToDeg();
            for (int b = 0; b < NB; ++b)
               if (thD >= blo[b] && thD < blo[b] + bw) {
                  double ex = dv_ex(m1, m2, m2, m1, Eb, thT, keN);
                  if (std::isfinite(ex))
                     h[b]->Fill(ex);
                  break;
               }
         }
         double mn = 1e9, mx = -1e9, sum = 0;
         int nOK = 0;
         printf("  %5.3f    %5.3f |", k, k * dvUsed);
         for (int b = 0; b < NB; ++b) {
            double pk = NAN;
            if (h[b]->Integral() > 120) {
               h[b]->Smooth(1);
               pk = h[b]->GetBinCenter(h[b]->GetMaximumBin());
            }
            printf(" %+5.2f", pk);
            if (std::isfinite(pk)) {
               mn = std::min(mn, pk);
               mx = std::max(mx, pk);
               sum += pk;
               ++nOK;
            }
            delete h[b];
         }
         double spread = (nOK >= 6) ? mx - mn : NAN;
         printf(" | %6.3f%s\n", spread, nOK < 6 ? "  (few bands)" : "");
         if (std::isfinite(spread)) {
            gSpread[i].SetPoint(npg[i], k * dvUsed, spread);
            gMuAvg[i].SetPoint(npg[i], k * dvUsed, sum / nOK);
            ++npg[i];
            if (spread < bestS[i]) {
               bestS[i] = spread;
               bestK[i] = k;
            }
         }
      }
      printf("  --> flattest at k = %.3f, i.e. dv = %.3f cm/us (spread %.3f MeV)\n", bestK[i], bestK[i] * dvUsed,
             bestS[i]);
   }

   TCanvas *c1 = new TCanvas("c1", "dv emulation", 1200, 500);
   c1->Divide(2, 1);
   for (int i = 0; i < NS; ++i) {
      gSpread[i].SetMarkerStyle(i == 0 ? 20 : 21);
      gSpread[i].SetMarkerColor(col[i]);
      gSpread[i].SetLineColor(col[i]);
      gSpread[i].SetLineWidth(2);
      gMuAvg[i].SetMarkerStyle(i == 0 ? 20 : 21);
      gMuAvg[i].SetMarkerColor(col[i]);
      gMuAvg[i].SetLineColor(col[i]);
      gMuAvg[i].SetLineWidth(2);
   }
   c1->cd(1);
   gSpread[0].SetTitle(TString::Format("E_{x} locus spread over #theta_{lab} 24-78 (E_{beam}=%.0f);"
                                       "emulated drift velocity [cm/#mus];max-min of the E_{x} peak [MeV]",
                                       Eb));
   gSpread[0].Draw("ALP");
   gSpread[1].Draw("LP same");
   auto *lg = new TLegend(0.62, 0.74, 0.89, 0.89);
   lg->AddEntry(&gSpread[0], lbl[0], "lp");
   lg->AddEntry(&gSpread[1], lbl[1], "lp");
   lg->Draw();
   c1->cd(2);
   gMuAvg[0].SetTitle("mean E_{x} of the elastic locus;emulated drift velocity [cm/#mus];<E_{x}> [MeV]");
   gMuAvg[0].Draw("ALP");
   gMuAvg[1].Draw("LP same");
   auto *z = new TLine(kLo * dvUsed, 0, kHi * dvUsed, 0);
   z->SetLineStyle(2);
   z->SetLineColor(kGray + 2);
   z->Draw();
   lg->Draw();

   TString png = here + TString::Format("/plots/dv_emulate_Eb%.0f.png", Eb);
   c1->SaveAs(png);
   printf("\nwrote %s\n\n", png.Data());
}
