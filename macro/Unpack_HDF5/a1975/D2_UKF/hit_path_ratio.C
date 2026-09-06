/// @file hit_path_ratio.C
/// @brief How far from "already ordered along the track" is a reco's RAW HIT ARRAY?
///
/// Metric: (path length walking the hit array as stored) / (path length of a nearest-neighbour
/// walk over the same hits). 1.0 means the array is already in track order; large means scrambled.
/// This is the quantity that decides whether AtPRA::SetTrackInitialParameters' azimuth unwrapping
/// -- which accumulates -Y crossings IN HIT-ARRAY ORDER -- is meaningful, and hence whether
/// AtPRA::OrderHitsAlongTrack changes anything for a given production.
///
///   root -b -q 'hit_path_ratio.C("/mnt/f/.../reco.root", 400, 30)'
#include <cmath>
#include <vector>

void hit_path_ratio(TString reco, int maxTracks = 400, int minHits = 30)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   auto f = TFile::Open(reco);
   if (!f || f->IsZombie()) { printf("cannot open %s\n", reco.Data()); return; }
   auto *t = (TTree *)f->Get("cbmsim");
   TClonesArray *pe = nullptr;
   t->SetBranchAddress("AtPatternEvent", &pe);

   std::vector<double> ratios, nh;
   for (Long64_t i = 0; i < t->GetEntries() && (int)ratios.size() < maxTracks; ++i) {
      t->GetEntry(i);
      auto *ev = (AtPatternEvent *)pe->At(0);
      if (!ev) continue;
      for (auto &trk : ev->GetTrackCand()) {
         if ((int)ratios.size() >= maxTracks) break;
         auto &hits = trk.GetHitArray();
         const int n = hits.size();
         if (n < minHits) continue;
         std::vector<ROOT::Math::XYZPoint> P;
         for (auto &h : hits) P.push_back(h->GetPosition());

         double stored = 0;
         for (int k = 1; k < n; ++k) stored += (P[k] - P[k - 1]).R();

         // greedy nearest-neighbour walk from the highest-Z hit (OrderHitsAlongTrack's convention)
         int seed = 0;
         for (int k = 1; k < n; ++k) if (P[k].Z() > P[seed].Z()) seed = k;
         std::vector<bool> used(n, false);
         used[seed] = true;
         int cur = seed;
         double nn = 0;
         for (int step = 1; step < n; ++step) {
            double best = 1e18; int bi = -1;
            for (int j = 0; j < n; ++j) {
               if (used[j]) continue;
               double d = (P[j] - P[cur]).R();
               if (d < best) { best = d; bi = j; }
            }
            if (bi < 0) break;
            nn += best; used[bi] = true; cur = bi;
         }
         if (nn > 1e-6) { ratios.push_back(stored / nn); nh.push_back(n); }
      }
   }
   if (ratios.empty()) { printf("no tracks\n"); return; }
   std::vector<double> s = ratios; std::sort(s.begin(), s.end());
   auto q = [&](double p) { return s[(int)(p * (s.size() - 1))]; };
   double mean = 0; for (double v : ratios) mean += v; mean /= ratios.size();
   int bad = 0; for (double v : ratios) if (v > 1.5) ++bad;
   double mnh = 0; for (double v : nh) mnh += v; mnh /= nh.size();
   printf("\n  %s\n", reco.Data());
   printf("  tracks %zu   mean hits/track %.0f\n", ratios.size(), mnh);
   printf("  stored-path / NN-path :  median %.2f   mean %.2f   p10 %.2f   p90 %.2f   max %.2f\n",
          q(0.5), mean, q(0.1), q(0.9), s.back());
   printf("  ratio > 1.5 (scrambled) : %d  (%.0f%%)\n", bad, 100.0 * bad / ratios.size());
}
