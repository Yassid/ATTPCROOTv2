/// @file check_vertex_beam_Be10.C
/// @brief Read the 10Be beam energy AT THE REACTION VERTEX back out of MC truth, and with it the
///        (t,p) kinematics -- before any reconstruction is involved.
///
/// WHY. Two separate things need checking before a campaign is worth running:
///
/// 1. THE BEAM ENERGY THAT WAS ACTUALLY TRANSPORTED. Every AT-TPC sim macro hands the ion mass in
///    amu to AtTPCIonGenerator, whose FairIon mass parameter is documented as GeV. If FairIon took
///    that literally the beam would be transported with the wrong mass and the energy implied by
///    pz would be wrong by ~7 %. This macro does not argue about it: it measures.
///
/// 2. THE CONSTANT Ebeam THE ANALYSIS SHOULD USE. acceptance_C14.C and ex_res_C14_hf.C invert the
///    two-body kinematics with ONE beam energy for every event, but the real beam has already lost
///    energy by the time it reaches the vertex, and the loss depends on where the vertex is. The
///    right constant is the MEAN over the reaction vertices, and it is obtained here by inverting
///    the kinematics the other way: for each truth event the level is known exactly (resEx), so
///    (theta_p, KE_p, Ex) determines E_beam. That is a measurement, not an estimate.
///
///   root -b -q 'check_vertex_beam_Be10.C("/mnt/f/Be10_tp/sims_b285/gs_s9001_sim.root", 0.0)'

#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TKey.h"
#include "TLegend.h"
#include "TLine.h"
#include "TMath.h"
#include "TProfile.h"
#include "TRandom.h"
#include "TString.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TFitResult.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <map>
#include <vector>
static double cvb_omega2(double x, double y, double z)
{
   return std::sqrt(x * x + y * y + z * z - 2 * x * y - 2 * y * z - 2 * x * z);
}

/// Ex(E_beam) for the measured (theta_lab, KE) of the ejectile -- the same expression the analysis
/// macros use. Solved for E_beam by bisection below, since it is monotonic over the range of use.
static double cvb_ex(double m1, double m2, double m3, double m4, double K_proj, double thetalab, double K_eject)
{
   double Et1 = K_proj + m1, Et3 = K_eject + m3;
   double s = m1 * m1 + m2 * m2 + 2 * m2 * Et1;
   double u = m2 * m2 + m3 * m3 - 2 * m2 * Et3;
   double arg = (std::cos(thetalab) * cvb_omega2(s, m1 * m1, m2 * m2) * cvb_omega2(u, m2 * m2, m3 * m3) -
                 (s - m1 * m1 - m2 * m2) * (m2 * m2 + m3 * m3 - u)) /
                   (2 * m2 * m2) + s + u - m2 * m2;
   if (arg <= 0)
      return std::numeric_limits<double>::quiet_NaN();
   return std::sqrt(arg) - m4;
}

void check_vertex_beam_Be10(TString f, Double_t resEx = 0.0, Double_t EbeamGen = 115.0,
                            Double_t mBeamAmu = 10.0135341, Double_t mTgtAmu = 3.0160493,
                            Double_t mEjAmu = 1.0078250, Double_t mResAmu = 12.0269221)
{
   gSystem->Load("libAtSimulationData.so");
   TFile *fi = TFile::Open(f);
   if (!fi || fi->IsZombie()) { printf("\033[1;31mcannot open %s\033[0m\n", f.Data()); return; }
   TTree *t = (TTree *)fi->Get("cbmsim");
   if (!t) { printf("\033[1;31mno cbmsim\033[0m\n"); return; }

   const double u = 931.49401;
   const double m1 = mBeamAmu * u, m2 = mTgtAmu * u, m3 = mEjAmu * u, m4 = mResAmu * u + resEx;

   TClonesArray *mc = nullptr;
   t->SetBranchAddress("MCTrack", &mc);

   TH1D *hE = new TH1D("hE", ";E_{beam} at the vertex [MeV];reactions", 200, EbeamGen - 20, EbeamGen + 5);
   TH1D *hZ = new TH1D("hZ", ";vertex z [mm];reactions", 100, 0, 1000);
   TH2D *hK = new TH2D("hK", ";#theta_{lab} [deg];KE_{p} [MeV]", 180, 0, 180, 200, 0, 40);
   TProfile *pEz = new TProfile("pEz", ";vertex z [mm];E_{beam} [MeV]", 50, 0, 1000);

   Long64_t N = t->GetEntries();
   long nreact = 0, nsolved = 0;
   double sumE = 0, sumE2 = 0;
   for (Long64_t i = 0; i < N; ++i) {
      t->GetEntry(i);
      if (!mc) continue;
      double keT = -1, thT = -1, zT = -1e9;
      for (int k = 0; k < mc->GetEntriesFast(); ++k) {
         auto *tr = (AtMCTrack *)mc->At(k);
         if (!tr || tr->GetPdgCode() != 2212 || tr->GetMotherId() != -1) continue;
         zT = tr->GetStartZ() * 10.0;
         double px = tr->GetPx() * 1000, py = tr->GetPy() * 1000, pz = tr->GetPz() * 1000;
         double p = std::sqrt(px * px + py * py + pz * pz);
         if (p <= 0) continue;
         keT = std::sqrt(p * p + m3 * m3) - m3;
         thT = std::acos(pz / p);
         break;
      }
      if (keT <= 0) continue;
      ++nreact;
      hZ->Fill(zT);
      hK->Fill(thT * TMath::RadToDeg(), keT);
      // bisect E_beam so that the reconstructed Ex equals the level that was generated
      double lo = EbeamGen - 40, hi = EbeamGen + 10;
      double flo = cvb_ex(m1, m2, m3, mResAmu * u, lo, thT, keT) - resEx;
      double fhi = cvb_ex(m1, m2, m3, mResAmu * u, hi, thT, keT) - resEx;
      if (std::isnan(flo) || std::isnan(fhi) || flo * fhi > 0) continue;
      for (int it = 0; it < 60; ++it) {
         double mid = 0.5 * (lo + hi);
         double fm = cvb_ex(m1, m2, m3, mResAmu * u, mid, thT, keT) - resEx;
         if (std::isnan(fm)) break;
         if (flo * fm <= 0) { hi = mid; fhi = fm; } else { lo = mid; flo = fm; }
      }
      double Ev = 0.5 * (lo + hi);
      ++nsolved;
      sumE += Ev; sumE2 += Ev * Ev;
      hE->Fill(Ev);
      pEz->Fill(zT, Ev);
   }
   double mean = nsolved ? sumE / nsolved : 0;
   double rms = nsolved ? std::sqrt(std::max(0.0, sumE2 / nsolved - mean * mean)) : 0;
   printf("\n================ 10Be(t,p)12Be truth check : %s ================\n", f.Data());
   printf("entries %lld, reactions %ld (%.1f %%), E_beam solved for %ld of them\n", N, nreact,
          100.0 * nreact / std::max<Long64_t>(1, N), nsolved);
   printf("level Ex = %.3f MeV, generated beam energy at the window = %.2f MeV\n", resEx, EbeamGen);
   printf("\033[1;32mE_beam AT THE VERTEX : mean %.3f MeV, rms %.3f, median-ish range [%.2f, %.2f]\033[0m\n", mean,
          rms, hE->GetXaxis()->GetBinLowEdge(hE->FindFirstBinAbove(0)),
          hE->GetXaxis()->GetBinUpEdge(hE->FindLastBinAbove(0)));
   printf("  -> pass this as Ebeam to acceptance_C14.C / ex_res_C14_hf.C\n");
   printf("  vertex z: mean %.1f mm, rms %.1f\n", hZ->GetMean(), hZ->GetRMS());
   printf("  loss from the window to the mean vertex: %.2f MeV\n", EbeamGen - mean);
   printf("\n  proton kinematics, KE vs theta_lab (median KE per band):\n");
   for (int b = 0; b < 12; ++b) {
      double a0 = b * 15.0, a1 = a0 + 15.0;
      TH1D *p = hK->ProjectionY("_py", hK->GetXaxis()->FindBin(a0 + 0.5), hK->GetXaxis()->FindBin(a1 - 0.5));
      if (p->GetEntries() < 5) { delete p; continue; }
      double q = 0.5, med;
      p->GetQuantiles(1, &med, &q);
      printf("    theta_lab %3.0f-%3.0f : n = %6.0f  median KE = %6.2f MeV\n", a0, a1, p->GetEntries(), med);
      delete p;
   }
   TCanvas *c = new TCanvas("c", "beam check", 1400, 500);
   c->Divide(3, 1);
   c->cd(1); hE->Draw();
   c->cd(2); pEz->Draw();
   c->cd(3); hK->Draw("colz");
   gSystem->mkdir("diagnostics", kTRUE);
   c->SaveAs("diagnostics/check_vertex_beam_Be10.png");
   printf("\nwrote diagnostics/check_vertex_beam_Be10.png\n");
   printf("beam check done\n");
}
