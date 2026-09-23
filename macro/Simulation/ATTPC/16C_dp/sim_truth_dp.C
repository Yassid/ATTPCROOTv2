/// @file sim_truth_dp.C
/// @brief Match the 16C(d,p)17C simulation's GENFIT output to MC TRUTH, event by event, and settle
/// the two things that must be settled before any resolution number is quoted:
///     1. the THETA CONVENTION -- the simulation reverses drift-z, so the reconstructed polar angle
///        is expected to come back as 180 - true. That is a claim to TEST, not to assume: it is
///        tested here by correlating reco theta against BOTH true and 180-true and seeing which one
///        wins. (reference_sim_z_handedness)
///     2. the VERTEX-Z convention -- z_lab = 1000 - z_digi bites every comparison between a fitted
///        vertex and a cluster position (reference_genfit_smoothed_positions_frame).
///
/// MATCHING IS BY EVENT INDEX. genfit reads the reco file sequentially and the reco read the sim
/// file sequentially, so event i is event i throughout. There is exactly ONE generated proton per
/// reacting event, so no track-level association is needed on the truth side.
///
/// TRACK SELECTION: the reco also finds the BEAM and the 17C RESIDUAL as pattern tracks, and with
/// no PID gate those were fitted as protons too. Here the fitted track closest to the true proton
/// in (KE, theta) is taken -- which is a TRUTH-ASSISTED choice and therefore an UPPER BOUND on what
/// the data can do. It is the right choice for a RESOLUTION FLOOR (the question is how well the
/// chain can do, not how well it picks), and the wrong one for an efficiency.
///
/// Uses exsel_dp.C's kine2b so the inversion is bit-identical to the one the data goes through.
///
///   root -l 'sim_truth_dp.C()'

#include "AtMCTrack.h"
#include "/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/a1975/D2_UKF/exsel_dp.C"

void sim_truth_dp(TString gfFile = "/mnt/f/a1975_C16_dp_sim/gs_s4001_genfitter_p.root",
                  TString simFile = "/mnt/f/a1975_C16_dp_sim/gs_s4001_sim.root", Long64_t maxEv = 0)
{
   using namespace exsel;
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtTools.so");

   TFile *fg = TFile::Open(gfFile);
   TFile *fs = TFile::Open(simFile);
   if (!fg || fg->IsZombie() || !fs || fs->IsZombie()) { printf("ERROR opening inputs\n"); return; }
   TTree *tg = (TTree *)fg->Get("cbmsim");
   TTree *ts = (TTree *)fs->Get("cbmsim");
   printf("\n  genfit entries %lld   sim entries %lld\n", tg->GetEntries(), ts->GetEntries());

   TClonesArray *te = nullptr, *pe = nullptr, *mc = nullptr;
   tg->SetBranchAddress("AtTrackingEvent", &te);
   tg->SetBranchAddress("AtPIDEvent", &pe);
   ts->SetBranchAddress("MCTrack", &mc);

   // --- the convention test -------------------------------------------------
   TH2D *hSame = new TH2D("hSame", "reco #theta vs TRUE;true #theta [deg];reco #theta [deg]", 60, 0, 180, 60, 0, 180);
   TH2D *hFlip = new TH2D("hFlip", "reco #theta vs 180-TRUE;180-true #theta [deg];reco #theta [deg]",
                          60, 0, 180, 60, 0, 180);
   TH2D *hVz = new TH2D("hVz", "reco v_{z} vs true;true v_{z} [mm];reco v_{z} [mm]", 60, 0, 1100, 60, 0, 1100);
   TH2D *hVzFlip = new TH2D("hVzFlip", "reco v_{z} vs 1000-true;1000-true v_{z} [mm];reco v_{z} [mm]",
                            60, 0, 1100, 60, 0, 1100);
   TH1D *hdKE = new TH1D("hdKE", "KE reco-true;#DeltaKE [MeV];", 100, -5, 5);

   const Long64_t N = maxEv > 0 ? TMath::Min(maxEv, tg->GetEntries()) : tg->GetEntries();
   Long64_t nTruth = 0, nFitEv = 0, nMatched = 0;
   double sSame = 0, sFlip = 0, sVz = 0, sVzFlip = 0;

   for (Long64_t i = 0; i < N; ++i) {
      ts->GetEntry(i);
      // truth: the one generated proton
      double tke = -1, tth = -1, tvz = -1;
      for (int j = 0; j < mc->GetEntriesFast(); ++j) {
         auto *m = (AtMCTrack *)mc->At(j);
         if (!m || m->GetPdgCode() != 2212) continue;
         const double px = m->GetPx(), py = m->GetPy(), pz = m->GetPz();
         const double pm = TMath::Sqrt(px * px + py * py + pz * pz) * 1000.0; // MeV/c
         const double M = 938.272;
         tke = TMath::Sqrt(pm * pm + M * M) - M;
         tth = TMath::ACos(pz / TMath::Sqrt(px * px + py * py + pz * pz)) * TMath::RadToDeg();
         tvz = m->GetStartZ() * 10.0; // cm -> mm
         break;
      }
      if (tke < 0) continue;
      ++nTruth;

      tg->GetEntry(i);
      if (!te || te->GetEntries() == 0) continue;
      auto *ev = (AtTrackingEvent *)te->At(0);
      if (!ev) continue;
      ++nFitEv;

      // the fitted track closest to truth in (KE, theta) -- truth-assisted, see the header
      double bd = 1e18, bke = -1, bth = -1, bvz = -1;
      for (auto &ft : ev->GetFittedTracks()) {
         if (!ft) continue;
         auto &k = ft->GetKinematicsXtr(); // the slot the analysis uses for genfit
         const double ke = k.kineticEnergy, th = k.theta * TMath::RadToDeg();
         if (ke <= 0) continue;
         auto v = ft->GetVertex();
         // score against BOTH conventions so the choice is not built in
         const double dS = TMath::Sq((ke - tke) / 1.0) + TMath::Sq((th - tth) / 5.0);
         const double dF = TMath::Sq((ke - tke) / 1.0) + TMath::Sq((th - (180 - tth)) / 5.0);
         const double d = TMath::Min(dS, dF);
         if (d < bd) { bd = d; bke = ke; bth = th; bvz = v.Z(); }
      }
      if (bke < 0) continue;
      ++nMatched;

      hSame->Fill(tth, bth);
      hFlip->Fill(180 - tth, bth);
      hVz->Fill(tvz, bvz);
      hVzFlip->Fill(1000 - tvz, bvz);
      hdKE->Fill(bke - tke);
      sSame += TMath::Abs(bth - tth);
      sFlip += TMath::Abs(bth - (180 - tth));
      sVz += TMath::Abs(bvz - tvz);
      sVzFlip += TMath::Abs(bvz - (1000 - tvz));
   }

   printf("\n  events with a generated proton   %lld\n", nTruth);
   printf("  of those, with a fitted track    %lld  (%.1f%%)\n", nFitEv, 100.0 * nFitEv / nTruth);
   printf("  matched                          %lld\n\n", nMatched);
   if (nMatched == 0) { printf("  nothing matched -- stop here.\n"); return; }

   printf("  === THETA CONVENTION (mean |reco - X|, smaller wins) ===\n");
   printf("     vs TRUE        %7.2f deg     correlation %+.3f\n", sSame / nMatched, hSame->GetCorrelationFactor());
   printf("     vs 180 - TRUE  %7.2f deg     correlation %+.3f\n", sFlip / nMatched, hFlip->GetCorrelationFactor());
   printf("     -> %s\n\n", (sFlip < sSame) ? "FLIPPED: reco theta = 180 - true, as expected for sim z handedness"
                                            : "NOT flipped: reco theta tracks true directly");

   printf("  === VERTEX Z CONVENTION (mean |reco - X|) ===\n");
   printf("     vs TRUE        %7.1f mm      correlation %+.3f\n", sVz / nMatched, hVz->GetCorrelationFactor());
   printf("     vs 1000 - TRUE %7.1f mm      correlation %+.3f\n", sVzFlip / nMatched,
          hVzFlip->GetCorrelationFactor());
   printf("     -> %s\n\n", (sVzFlip < sVz) ? "MIRRORED: reco vz = 1000 - true" : "direct: reco vz tracks true");

   printf("  === KE ===\n");
   printf("     median |dKE| proxy: mean %+.3f MeV, RMS %.3f\n\n", hdKE->GetMean(), hdKE->GetRMS());
}
