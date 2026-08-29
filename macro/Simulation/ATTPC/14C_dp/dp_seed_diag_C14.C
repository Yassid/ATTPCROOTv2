/// @file dp_seed_diag_C14.C
/// @brief Is the backward-track problem in the SEED or in the FIT?
///
///   root -b -q 'dp_seed_diag_C14.C()'
///
/// AtGenfitter builds its seed from two quantities the pattern recognition supplies:
///     GeoRadius  -- the radius of the circle fitted to the hits projected on the pad plane
///     GeoTheta   -- the polar angle, from a RANSAC line in the (arc length, z) plane, and the
///                   ONLY thing deciding whether a track is treated as backward
/// and then p = 0.3 * |B| * R / sin(theta). Everything downstream inherits that. This macro puts
/// the seed, the fit and the truth side by side on the same tracks, so the failure can be
/// attributed instead of guessed.
///
/// GeoTheta deserves particular suspicion: AtPatternBridgeTask computes it as acos(|dz|/r3), which
/// cannot exceed 90 deg at all, while AtPRA computes acos(sign*|dir.Y|) which can. Only the second
/// is in this chain, but the hemisphere it returns rests on the sign of a RANSAC direction in a
/// plane where a backward track and a forward one differ only by that sign.

#include <algorithm>
#include <vector>

static double sd_q(std::vector<double> v, double p)
{
   if (v.size() < 10) return NAN;
   size_t k = (size_t)std::min<double>(v.size() - 1, std::max(0.0, p * (v.size() - 1)));
   std::nth_element(v.begin(), v.begin() + k, v.end());
   return v[k];
}

void dp_seed_diag_C14(TString simFile = "/mnt/f/a1954_C14dp_hf/sims_b285/gs_s8001_sim.root",
                      TString recoFile = "/mnt/f/a1954_C14dp_matfx/ab_reco.root",
                      TString fitFile = "/mnt/f/a1954_C14dp_matfx/matON_genfit.root", Double_t bField = 2.85)
{
   gSystem->Load("libAtReconstruction.so");
   gSystem->Load("libAtSimulationData.so");
   const double u = 931.49401, mp = 1.007825 * u;

   TFile *fs = TFile::Open(simFile);
   TTree *ts = (TTree *)fs->Get("cbmsim");
   TClonesArray *mc = nullptr;
   ts->SetBranchAddress("MCTrack", &mc);

   TFile *fr = TFile::Open(recoFile);
   TTree *tr = (TTree *)fr->Get("cbmsim");
   TClonesArray *pe = nullptr;
   tr->SetBranchAddress("AtPatternEvent", &pe);

   TFile *ff = TFile::Open(fitFile);
   TTree *tf = ff ? (TTree *)ff->Get("cbmsim") : nullptr;
   TClonesArray *te = nullptr;
   if (tf) tf->SetBranchAddress("AtTrackingEvent", &te);

   const int NB = 5;
   const double lo[NB] = {30, 60, 85, 95, 115}, hi[NB] = {60, 85, 95, 115, 150};
   std::vector<double> dTheta[NB], rRatio[NB], pSeed[NB], pFit[NB], nHit[NB];
   long nBackTagged[NB] = {0}, nTot[NB] = {0};

   Long64_t N = std::min(ts->GetEntries(), tr->GetEntries());
   for (Long64_t i = 0; i < N; ++i) {
      ts->GetEntry(i);
      if (!mc) continue;
      double keT = -1, thT = -1;
      for (int k = 0; k < mc->GetEntriesFast(); ++k) {
         auto *p = (AtMCTrack *)mc->At(k);
         if (!p || p->GetPdgCode() != 2212 || p->GetMotherId() != -1) continue;
         double px = p->GetPx() * 1000, py = p->GetPy() * 1000, pz = p->GetPz() * 1000;
         double pm = std::sqrt(px * px + py * py + pz * pz);
         if (pm <= 0) continue;
         keT = std::sqrt(pm * pm + mp * mp) - mp;
         thT = std::acos(pz / pm) * TMath::RadToDeg();
         break;
      }
      if (keT <= 0) continue;
      int b = -1;
      for (int j = 0; j < NB; ++j) if (thT >= lo[j] && thT < hi[j]) b = j;
      if (b < 0) continue;
      ++nTot[b];

      tr->GetEntry(i);
      if (!pe || pe->GetEntriesFast() == 0) continue;
      auto *ev = (AtPatternEvent *)pe->At(0);
      if (!ev) continue;
      auto &cand = ev->GetTrackCand();
      if (cand.empty()) continue;
      // the longest track in the event -- with the beam hole inhibited there is essentially one
      const AtTrack *best = nullptr;
      size_t nb = 0;
      for (auto &t : cand)
         if (t.GetHitArray().size() > nb) { nb = t.GetHitArray().size(); best = &t; }
      if (!best) continue;

      double geoTh = best->GetGeoTheta() * TMath::RadToDeg();
      double geoR = best->GetGeoRadius(); // mm
      if (!(geoR > 0)) continue;
      // truth transverse radius for this proton
      double pTrue = std::sqrt(keT * keT + 2 * keT * mp) / 1000.0;              // GeV/c
      double rTrue = 1000.0 * pTrue * std::sin(thT * TMath::DegToRad()) / (0.299792458 * bField); // mm
      // the seed momentum genfit would build from this pattern, using GeoTheta as the polar
      double sth = std::abs(std::sin(geoTh * TMath::DegToRad()));
      if (sth < 1e-3) sth = 1e-3;
      double brho = bField * (geoR / 1000.0) / sth;
      double pS = 0.299792458 * brho * 1000.0; // MeV/c
      double keS = std::sqrt(pS * pS + mp * mp) - mp;

      dTheta[b].push_back(geoTh - thT);
      rRatio[b].push_back(geoR / rTrue);
      pSeed[b].push_back(keS / keT);
      nHit[b].push_back(nb);
      if (geoTh > 90) ++nBackTagged[b];

      if (tf) {
         tf->GetEntry(i);
         if (te && te->GetEntriesFast() > 0) {
            auto *fe = (AtTrackingEvent *)te->At(0);
            if (fe)
               for (auto &ft : fe->GetFittedTracks()) {
                  if (!ft) continue;
                  const auto &md = ft->GetTrackMetadata();
                  double ndf = md ? md->GetNdf() : 0, c2 = md ? md->GetChi2() : 0;
                  if (!(ndf > 0 && c2 / ndf < 5)) continue;
                  double ke = ft->GetKinematicsXtr().kineticEnergy;
                  if (ke > 0 && ke < 1000) pFit[b].push_back(ke / keT);
                  break;
               }
         }
      }
   }

   printf("\n  14C(d,p) at %.2f T -- seed quality against truth, by TRUE proton angle\n", bField);
   printf("  %-10s %6s %7s | %10s %9s | %10s %9s | %9s %9s | %s\n", "theta_lab", "n", "<hits>", "GeoTheta-true",
          "IQR", "R_geo/R_true", "IQR", "KE_seed/true", "KE_fit/true", "tagged backward");
   for (int b = 0; b < NB; ++b) {
      if (dTheta[b].size() < 10) { printf("  %3.0f-%-6.0f %6zu   (too few)\n", lo[b], hi[b], dTheta[b].size()); continue; }
      printf("  %3.0f-%-6.0f %6zu %7.0f | %10.1f %9.1f | %10.3f %9.3f | %9.3f %9.3f | %5.0f %%\n", lo[b], hi[b],
             dTheta[b].size(), sd_q(nHit[b], .5), sd_q(dTheta[b], .5), (sd_q(dTheta[b], .75) - sd_q(dTheta[b], .25)),
             sd_q(rRatio[b], .5), (sd_q(rRatio[b], .75) - sd_q(rRatio[b], .25)), sd_q(pSeed[b], .5),
             pFit[b].size() >= 10 ? sd_q(pFit[b], .5) : NAN, 100.0 * nBackTagged[b] / dTheta[b].size());
   }
   printf("\n  R_geo/R_true is the pattern circle against the true transverse radius: the seed momentum\n"
          "  is proportional to it. KE_seed/true is what that radius plus GeoTheta implies before genfit\n"
          "  touches it; KE_fit/true is what comes out. 'tagged backward' is GeoTheta > 90 deg, the only\n"
          "  thing that turns the backward seed on.\n\n");
   fs->Close(); fr->Close(); if (ff) ff->Close();
}
