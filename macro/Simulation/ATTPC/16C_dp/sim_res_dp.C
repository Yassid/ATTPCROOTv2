/// @file sim_res_dp.C
/// @brief The RESOLUTION FLOOR of 16C(d,p)17C, vs VERTEX Z -- the measurement the simulation was
/// built for.
///
/// THE QUESTION. Under the adopted data selection the BOUND-region low-Ex peak walks +0.30 MeV
/// across the chamber and its width collapses 4x at v_z 650-750 mm
/// ([[project_a1975_dp_vz_drift]]). Both are either geometry (the detector really does this, and
/// the analysis must correct for it) or reconstruction (a defect that could be fixed). Truth
/// through the same chain is the only thing that separates them.
///
/// WHAT IS MEASURED, per vertex-z slice:
///     dKE   = KE_reco  - KE_true            <- is the fitter's energy bias vz-dependent?
///     dEx   = Ex(reco) - Ex(truth)          <- the SAME inversion both sides, so the level
///                                              structure and the beam energy cancel exactly
///     sigma(dEx)                            <- the resolution FLOOR
/// MEDIANS and quantile widths, never means: a fit-failure tail makes a mean say anything
/// ([[feedback_use_medians]]). sigma is quoted as the half 16-84 % interval, which is the gaussian
/// sigma for a gaussian and is not hijacked by outliers.
///
/// NO keOff IS APPLIED. The data selection carries an empirical +0.25 MeV KE offset; the whole
/// point here is to measure what the offset should be, so applying it would build in the answer.
/// The quality cuts (theta window, chi2, vz) ARE applied, so the sample is the analysed one.
///
/// TRUTH-ASSISTED TRACK CHOICE, and it is an UPPER BOUND. With no PID gate the beam and the 17C
/// residual were fitted as protons too; the track closest to truth is taken. That measures how
/// well the chain CAN do, which is what a floor is. It is NOT an efficiency and must not be read
/// as one.
///
///   root -l 'sim_res_dp.C()'

#include "AtMCTrack.h"
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"

namespace {
double vmed(std::vector<double> v)
{
   if (v.empty()) return -999;
   std::sort(v.begin(), v.end());
   const size_t n = v.size();
   return (n % 2) ? v[n / 2] : 0.5 * (v[n / 2 - 1] + v[n / 2]);
}
double vq(std::vector<double> v, double q)
{
   if (v.empty()) return -999;
   std::sort(v.begin(), v.end());
   return v[std::min(v.size() - 1, (size_t)(q * v.size()))];
}
/// half the 16-84 % interval: the gaussian sigma for a gaussian, robust against tails
double vsig(std::vector<double> v)
{
   if (v.size() < 5) return -999;
   return 0.5 * (vq(v, 0.84) - vq(v, 0.16));
}
} // namespace

void sim_res_dp(TString gfFile = "/mnt/f/a1975_C16_dp_sim/gs_s4001_genfitter_p.root",
                TString simFile = "/mnt/f/a1975_C16_dp_sim/gs_s4001_sim.root", Long64_t maxEv = 0,
                Double_t chi2Cut = 13.5, Double_t thLo = 89, Double_t thHi = 161, Double_t keLo = 1.2)
{
   using namespace exsel;
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");
   Cuts s;
   s.keOff = 0.0; // DELIBERATE -- see the header
   s.kcOn = false;

   TFile *fg = TFile::Open(gfFile);
   TFile *fs = TFile::Open(simFile);
   TTree *tg = (TTree *)fg->Get("cbmsim");
   TTree *ts = (TTree *)fs->Get("cbmsim");
   TClonesArray *te = nullptr, *mc = nullptr;
   tg->SetBranchAddress("AtTrackingEvent", &te);
   ts->SetBranchAddress("MCTrack", &mc);

   const int NS = 8;
   const double zlo = 150, zw = 100;
   std::vector<double> dKE[NS], dEx[NS], exT[NS], allDKE, allDEx;
   long nsl[NS] = {0};
   Long64_t nTruth = 0, nPass = 0;

   const Long64_t N = maxEv > 0 ? TMath::Min(maxEv, tg->GetEntries()) : tg->GetEntries();
   for (Long64_t i = 0; i < N; ++i) {
      ts->GetEntry(i);
      double tke = -1, tth = -1, tvz = -1;
      for (int j = 0; j < mc->GetEntriesFast(); ++j) {
         auto *m = (AtMCTrack *)mc->At(j);
         if (!m || m->GetPdgCode() != 2212) continue;
         const double px = m->GetPx(), py = m->GetPy(), pz = m->GetPz();
         const double pm = TMath::Sqrt(px * px + py * py + pz * pz) * 1000.0;
         const double M = 938.272;
         tke = TMath::Sqrt(pm * pm + M * M) - M;
         tth = TMath::ACos(pz / TMath::Sqrt(px * px + py * py + pz * pz)) * TMath::RadToDeg();
         tvz = m->GetStartZ() * 10.0;
         break;
      }
      if (tke < 0) continue;
      ++nTruth;

      tg->GetEntry(i);
      if (!te || te->GetEntries() == 0) continue;
      auto *ev = (AtTrackingEvent *)te->At(0);
      if (!ev) continue;

      double bd = 1e18, bke = -1, bth = -1, bvz = -1, bc2 = 1e9;
      for (auto &ft : ev->GetFittedTracks()) {
         if (!ft) continue;
         auto &k = ft->GetKinematicsXtr();
         const double ke = k.kineticEnergy, th = k.theta * TMath::RadToDeg();
         if (ke <= 0) continue;
         // GetTrackMetadata returns a unique_ptr REFERENCE, not a raw pointer -- `auto *` does not
         // bind to it. ex_dp_dv1104.C calls through it directly for the same reason.
         auto &md = ft->GetTrackMetadata();
         const double ndf = md ? md->GetNdf() : -1, c2 = md ? md->GetChi2() : 0;
         const double c2n = ndf > 0 ? c2 / ndf : 1e9;
         auto v = ft->GetVertex();
         const double d = TMath::Sq((ke - tke) / 1.0) + TMath::Sq((th - tth) / 5.0);
         if (d < bd) { bd = d; bke = ke; bth = th; bvz = v.Z(); bc2 = c2n; }
      }
      if (bke < 0) continue;

      // the ANALYSIS cuts, on the reconstructed quantities, exactly as the data is cut
      if (!(bc2 <= chi2Cut)) continue;
      if (!(bth >= thLo && bth <= thHi)) continue;
      if (!(bke >= keLo)) continue;
      if (!(bvz >= -100 && bvz <= 940)) continue;
      ++nPass;

      double exR = 0, tcR = 0, exT_ = 0, tcT = 0;
      if (!kine2b(s.ebeam, bth * TMath::Pi() / 180.0, bke, exR, tcR)) continue;
      if (!kine2b(s.ebeam, tth * TMath::Pi() / 180.0, tke, exT_, tcT)) continue;

      const int b = (int)((bvz - zlo) / zw);
      allDKE.push_back(bke - tke);
      allDEx.push_back(exR - exT_);
      if (b >= 0 && b < NS) {
         dKE[b].push_back(bke - tke);
         dEx[b].push_back(exR - exT_);
         exT[b].push_back(exT_);
         ++nsl[b];
      }
   }

   printf("\n=== sim_res_dp : the RESOLUTION FLOOR, 16C(d,p)17C g.s. ===\n");
   printf("  cuts: chi2/ndf <= %.1f, theta %.0f-%.0f, KE >= %.1f, vz -100..940.  NO keOff.\n",
          chi2Cut, thLo, thHi, keLo);
   printf("  events with a generated proton %lld   passing the analysis cuts %lld (%.1f%%)\n\n",
          nTruth, nPass, 100.0 * nPass / TMath::Max((Long64_t)1, nTruth));
   if (nPass < 20) { printf("  too few to say anything.\n"); return; }

   printf("  ALL SLICES TOGETHER\n");
   printf("    median dKE  = %+.3f MeV   (sigma %.3f)   <- compare the analysis keOff = +0.25\n",
          vmed(allDKE), vsig(allDKE));
   printf("    median dEx  = %+.3f MeV   (sigma %.3f)   <- THE RESOLUTION FLOOR\n\n",
          vmed(allDEx), vsig(allDEx));

   printf("  BY VERTEX Z -- the question: does the FLOOR drift, and does it collapse at 650-750?\n");
   printf("    vz slice      n    median dKE   median dEx   sigma(dEx)   true Ex\n");
   for (int b = 0; b < NS; ++b) {
      if (nsl[b] < 10) { printf("    %4.0f-%4.0f   %4ld    (too few)\n", zlo + zw * b, zlo + zw * (b + 1), nsl[b]); continue; }
      printf("    %4.0f-%4.0f   %4ld     %+7.3f      %+7.3f      %6.3f     %+6.3f\n", zlo + zw * b,
             zlo + zw * (b + 1), nsl[b], vmed(dKE[b]), vmed(dEx[b]), vsig(dEx[b]), vmed(exT[b]));
   }
   printf("\n  DATA, for comparison (project_a1975_dp_vz_drift, peak Ex per slice):\n");
   printf("    150-250 -0.190 | 250-350 -0.127 | 350-450 -0.119 | 450-550 -0.070\n");
   printf("    550-650 +0.052 | 650-750 -0.020 (sigma 0.927!) | 750-850 +0.110 | 850-950 +0.055\n");
   printf("  The data walks +0.30 MeV end to end. If this simulation walks the same way, the drift\n");
   printf("  is the DETECTOR and must be corrected; if it is flat, the drift is in the DATA and is\n");
   printf("  something the simulation does not contain (calibration, beam energy, gas, or dv).\n\n");
}
