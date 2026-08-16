/// @file beam_axis_check.C
/// @brief Is the beam actually on the geometric z-axis the back-extrapolation assumes?
///
/// AtGenfitter extrapolates the fitted state to the LINE (0,0,0)+(0,0,1), so the stored vertex
/// is the point of closest approach to the ASSUMED axis and its (x,y) is the impact-parameter
/// vector. Nothing measures the beam position, but the tracks constrain it:
///
///   - beam ON the assumed axis: tracks originate on it, so the impact parameter is ~0 and its
///     (x,y) centroid is 0, with a spread that is just vertex resolution.
///   - beam OFFSET by d: backward-extrapolated tracks do not pass through the assumed axis. The
///     centroid moves to ~d, and vr picks up a MODULATION with the track azimuth, because a
///     track emitted towards the axis crosses it while one emitted away does not.
///
/// The modulation is the discriminator. A centroid alone can be moved by any left-right
/// asymmetry in acceptance; a phi-dependence that follows cos(phi - phi0) is the signature of a
/// transverse displacement specifically.
///
/// Also counts multi-track events, since a two-track closest approach is a vertex that assumes
/// NOTHING about the beam and is the only fully independent check available here.
///
///   root -b -q 'beam_axis_check.C("/mnt/f/a1975/gf_dt_ab_nodouble/run_0031_multifit_genfitter_t.root")'

#include <algorithm>
#include <cmath>
#include <vector>

void beam_axis_check(TString file)
{
   gSystem->Load("libAtReconstruction.so");
   TFile *f = TFile::Open(file);
   TTree *t = (TTree *)f->Get("cbmsim");
   TClonesArray *te = nullptr;
   t->SetBranchAddress("AtTrackingEvent", &te);

   std::vector<double> vx, vy, vr, vz;
   // vr in 8 bins of track azimuth: a transverse offset shows up as a cos(phi-phi0) modulation
   std::vector<std::vector<double>> byPhi(8);
   long nEvt = 0, nMulti = 0, nTrk = 0;

   for (Long64_t i = 0; i < t->GetEntries(); ++i) {
      t->GetEntry(i);
      if (!te || te->GetEntries() == 0)
         continue;
      auto *ev = (AtTrackingEvent *)te->At(0);
      if (!ev)
         continue;
      int nGood = 0;
      for (auto &ft : ev->GetFittedTracks()) {
         if (!ft)
            continue;
         const auto &m = ft->GetTrackMetadata();
         if (!m || !(m->GetNdf() > 0 && m->GetChi2() > 0))
            continue; // converged fits only
         ++nGood;
         auto v = ft->GetVertex(); // mm
         const double x = v.X(), y = v.Y(), z = v.Z();
         const double r = std::sqrt(x * x + y * y);
         if (!std::isfinite(r) || r > 200)
            continue; // runaway extrapolation
         vx.push_back(x);
         vy.push_back(y);
         vr.push_back(r);
         vz.push_back(z);
         double phi = ft->GetKinematicsXtr().phi; // rad
         if (!std::isfinite(phi))
            continue;
         while (phi < 0)
            phi += 2 * TMath::Pi();
         int b = std::min(7, (int)(8 * phi / (2 * TMath::Pi())));
         byPhi[b].push_back(r);
      }
      nTrk += nGood;
      if (nGood > 0)
         ++nEvt;
      if (nGood >= 2)
         ++nMulti;
   }

   auto med = [](std::vector<double> v) {
      if (v.empty())
         return 0.0;
      std::sort(v.begin(), v.end());
      return v[v.size() / 2];
   };
   auto qspread = [](std::vector<double> v) {
      if (v.size() < 5)
         return 0.0;
      std::sort(v.begin(), v.end());
      return 0.5 * (v[(size_t)(0.84 * v.size())] - v[(size_t)(0.16 * v.size())]);
   };

   printf("\n=== reconstructed vertex, converged fits only (n=%zu) ===\n", vx.size());
   printf("  median x = %+8.3f mm   (16-84 halfwidth %6.3f)\n", med(vx), qspread(vx));
   printf("  median y = %+8.3f mm   (16-84 halfwidth %6.3f)\n", med(vy), qspread(vy));
   printf("  median z = %+8.1f mm   (16-84 halfwidth %6.1f)\n", med(vz), qspread(vz));
   printf("  median r = %8.3f mm   (16-84 halfwidth %6.3f)\n", med(vr), qspread(vr));
   printf("\n  If the beam sat on the assumed axis, median r is vertex RESOLUTION.\n");
   printf("  A median r well above that, with the phi structure below, means an OFFSET.\n");

   printf("\n=== impact parameter vs track azimuth (offset -> cos(phi-phi0) modulation) ===\n");
   double lo = 1e9, hi = -1e9;
   for (int b = 0; b < 8; ++b) {
      const double m = med(byPhi[b]);
      printf("  phi %3d-%3d deg : n=%-5zu  median r = %7.3f mm\n", b * 45, (b + 1) * 45, byPhi[b].size(), m);
      if (byPhi[b].size() > 20) {
         lo = std::min(lo, m);
         hi = std::max(hi, m);
      }
   }
   if (hi > lo)
      printf("  modulation amplitude (max-min)/2 = %.3f mm\n", 0.5 * (hi - lo));

   printf("\n=== multi-track events (the axis-INDEPENDENT check) ===\n");
   printf("  events with >=1 converged track : %ld\n", nEvt);
   printf("  events with >=2 converged tracks: %ld  (%.1f%%)\n", nMulti, nEvt ? 100.0 * nMulti / nEvt : 0.0);
   printf("  -> a track-track closest approach needs >=2; %s\n\n",
          nMulti > 100 ? "enough here to try it" : "NOT enough in this run, would need the full set");
}
