/// @file pt_truth_check.C
/// @brief Validate the 16C(p,t)14C generation against the analytic locus, BEFORE anything is
/// digitised. If the truth is wrong, every downstream number is wrong and looks plausible.
///
/// Three things are checked, each of which has burned this analysis before:
///
///  1. THE LOCUS. The (p,t) locus is double-valued: theta_lab rises to about 37 deg near residual
///     cm 45-60 and folds BACK to 0.8 deg, so every lab angle carries a low-KE and a high-KE
///     solution. Unlike (d,t) the high branch is not obviously unmeasurable here (its cyclotron
///     radius peaks near 298 mm against a ~290 mm chamber, not 840), so BOTH are generated and the
///     median below is a mixture. Analytic LOW branch at 185 MeV: 5.0 MeV at 10 deg, 5.9 at 20,
///     8.2 at 30. The medians should sit at or above those.
///
///  2. THE VERTEX. maxELoss controls where the reaction fires. If it is wrong the vertex piles up
///     at the entrance and the acceptance is then measured on the wrong target thickness.
///
///  3. THE BEAM ENERGY AT THE VERTEX. Solving the two-body kinematics with the ENTRANCE energy
///     while the reaction really happened deeper in the gas puts the apparent Ex at the wrong
///     place, and the error grows with z. This is not small: on (p,d) it reached 1.3 MeV across
///     the chamber, twice the 0.740 MeV level spacing.
///
///   root -b -q 'pt_truth_check.C("./data/dt_test.root",0.0)'

#include <tuple>

static double ptc_omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

/// two-body kinematics: returns {Ex, theta_cm[deg]} for a measured ejectile (KE, theta_lab)
static std::tuple<double, double> ptc_kine(double m1, double m2, double m3, double m4, double K_proj, double thetalab,
                                           double K_eject)
{
   double Et1 = K_proj + m1, Et3 = K_eject + m3, Et4 = Et1 + m2 - Et3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double m4_ex = std::sqrt((std::cos(thetalab) * ptc_omega2(s, m1 * m1, m2 * m2) * ptc_omega2(u, m2 * m2, m3 * m3) -
                             (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                               (2 * m2 * m2) +
                            s + u - m2 * m2);
   double Ex = m4_ex - m4;
   double t = m2 * m2 + m4_ex * m4_ex - 2 * m2 * Et4;
   double theta_cm = TMath::Pi() - std::acos((s * s + s * (2 * t - m1 * m1 - m2 * m2 - m3 * m3 - m4_ex * m4_ex) +
                                              (m1 * m1 - m2 * m2) * (m3 * m3 - m4_ex * m4_ex)) /
                                             (ptc_omega2(s, m1 * m1, m2 * m2) * ptc_omega2(s, m3 * m3, m4_ex * m4_ex)));
   return {Ex, theta_cm * TMath::RadToDeg()};
}

void pt_truth_check(TString simFile = "./data/smoke_gs.root", Double_t resEx = 0.0, Double_t Ebeam = 185.0)
{
   gSystem->Load("libAtSimulationData.so");
   gStyle->SetOptStat(0);

   const double u = 931.49401;
   const double m_C16 = 16.0147 * u, m_p = 1.0078250322 * u;
   const double m_t = 3.0160492779 * u, m_C14 = 14.0032420 * u;
   const double m_resid = m_C14 + resEx;

   TFile *fs = TFile::Open(simFile);
   TTree *ts = fs && !fs->IsZombie() ? (TTree *)fs->Get("cbmsim") : nullptr;
   if (!ts) { printf("\n  cannot open %s\n\n", simFile.Data()); return; }
   TClonesArray *mc = nullptr;
   ts->SetBranchAddress("MCTrack", &mc);

   const int NL = 12; // 5 deg bins in theta_lab, 0-60
   std::vector<double> keInBin[NL];
   std::vector<double> exInZ[10]; // Ex vs vertex z, 100 mm slabs
   long nEv = 0, nReact = 0;
   double zMin = 1e9, zMax = -1e9;
   auto *hZ = new TH1D("hZ", "reaction vertex;z [mm];reactions", 40, 0, 1000);
   auto *hCM = new TH1D("hCM", "generated;#theta_{cm} [deg];reactions", 36, 0, 180);
   auto *hKEth = new TH2D("hKEth", "truth locus;#theta_{lab} [deg];KE_{t} [MeV]", 90, 0, 90, 100, 0, 90);

   const Long64_t N = ts->GetEntries();
   for (Long64_t i = 0; i < N; ++i) {
      ts->GetEntry(i);
      ++nEv;
      if (!mc) continue;
      double ke = -1, th = -1, z = -1e9;
      for (int k = 0; k < mc->GetEntriesFast(); ++k) {
         auto *tr = (AtMCTrack *)mc->At(k);
         if (!tr || tr->GetPdgCode() != 1000010030 || tr->GetMotherId() != -1) // triton, primary
            continue;
         z = tr->GetStartZ() * 10.0; // cm -> mm
         double px = tr->GetPx() * 1000, py = tr->GetPy() * 1000, pz = tr->GetPz() * 1000;
         double p = std::sqrt(px * px + py * py + pz * pz);
         if (p <= 0) continue;
         ke = std::sqrt(p * p + m_t * m_t) - m_t;
         th = std::acos(pz / p);
         break;
      }
      if (ke <= 0) continue; // beam-only event
      ++nReact;
      const double thDeg = th * TMath::RadToDeg();
      hZ->Fill(z);
      hKEth->Fill(thDeg, ke);
      zMin = std::min(zMin, z);
      zMax = std::max(zMax, z);
      int lb = (int)(thDeg / 5.0);
      if (lb >= 0 && lb < NL) keInBin[lb].push_back(ke);
      auto [ex, cm] = ptc_kine(m_C16, m_p, m_t, m_resid, Ebeam, th, ke);
      if (!std::isnan(cm)) hCM->Fill(cm);
      int zb = (int)(z / 100.0);
      if (zb >= 0 && zb < 10 && !std::isnan(ex)) exInZ[zb].push_back(ex);
   }

   auto med = [](std::vector<double> &v) {
      if (v.empty()) return std::numeric_limits<double>::quiet_NaN();
      std::sort(v.begin(), v.end());
      return v[v.size() / 2];
   };

   printf("\n  %lld events, %ld with a triton (%.1f %% reacted)\n", N, nReact, 100.0 * nReact / std::max(1L, nEv));
   printf("\n  ===== 1. THE LOCUS: median KE vs theta_lab, against the analytic low branch =====\n");
   printf("  %-12s %8s %10s   %s\n", "theta_lab", "N", "KE median", "expected (low branch)");
   // The (d,t) numbers do not transfer. Computed analytically for 16C(p,t)14C at 185 MeV on the
   // LOW-KE branch (residual cm angle below the ~60 deg turnover):
   const double expTh[3] = {10, 20, 30}, expKE[3] = {5.0, 5.9, 8.2};
   for (int b = 0; b < NL; ++b) {
      if (keInBin[b].empty()) continue;
      const double c = 5.0 * b + 2.5;
      TString note;
      for (int e = 0; e < 3; ++e)
         if (std::fabs(c - expTh[e]) < 2.6) note = Form("  <-- expect %.2f MeV", expKE[e]);
      printf("  %4.0f - %-4.0f %8zu %10.2f%s\n", 5.0 * b, 5.0 * (b + 1), keInBin[b].size(), med(keInBin[b]),
             note.Data());
   }

   printf("\n  ===== 2. THE VERTEX: must be flat along the drift length =====\n");
   printf("  z range %.0f to %.0f mm\n", zMin, zMax);
   double lo = hZ->Integral(1, hZ->GetNbinsX() / 2), hi = hZ->Integral(hZ->GetNbinsX() / 2 + 1, hZ->GetNbinsX());
   printf("  first half %.0f, second half %.0f, ratio %.2f  (1.00 = flat; a ratio far above 1\n"
          "  means maxELoss is too small and the reactions pile up at the entrance)\n",
          lo, hi, hi > 0 ? lo / hi : -1);

   printf("\n  ===== 3. BEAM ENERGY AT THE VERTEX: apparent Ex vs z, solved with the ENTRANCE energy =====\n");
   printf("  (the residual was generated at Ex = %.3f MeV; drift is the beam energy loss)\n", resEx);
   printf("  %-14s %8s %10s %10s\n", "z slab [mm]", "N", "Ex median", "drift");
   double ex0 = std::numeric_limits<double>::quiet_NaN();
   for (int zb = 0; zb < 10; ++zb) {
      if (exInZ[zb].size() < 5) continue;
      double m = med(exInZ[zb]);
      if (std::isnan(ex0)) ex0 = m;
      printf("  %4d - %-7d %8zu %10.3f %+10.3f\n", 100 * zb, 100 * (zb + 1), exInZ[zb].size(), m, m - ex0);
   }

   auto *c1 = new TCanvas("c_truth", "truth check", 1400, 500);
   c1->Divide(3, 1);
   c1->cd(1); hKEth->Draw("colz");
   c1->cd(2); hZ->Draw("hist");
   c1->cd(3); hCM->Draw("hist");
   c1->SaveAs("./data/pt_truth_check.png");
   printf("\n  wrote data/pt_truth_check.png\n\n");
}
