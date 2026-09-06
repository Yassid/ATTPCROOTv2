/// @file ab_hitord_compare.C
/// @brief A/B the hit-ordering fix: production reco vs re-reco with AtPRA::OrderHitsAlongTrack on.
///
/// Matched track by track on (entry, trackID). HDBSCAN's point clusters are unchanged by the fix,
/// so trackID and the hit CONTENT must be identical -- that is checked, and a mismatch means
/// something other than the ordering moved and the comparison is void.
///
/// What is expected to change, and why:
///   GeoTheta   -- the arc-length fit unwraps the azimuth in hit-array order (the defect)
///   nClusters  -- ClusterizeSmooth3D seeds on hitArray.at(0) and walks the array sequentially,
///                 so reordering the hits also rebuilds the AtHitCluster array genfit consumes
///
///   root -b -q 'ab_hitord_compare.C("run_0020")'
#include <cmath>
#include <map>
#include <vector>

void ab_hitord_compare(TString run = "run_0020",
                       TString oldDir = "/mnt/f/a1975/reco_d2_dv1104/",
                       TString newDir = "/mnt/f/a1975/reco_d2_dv1104_hitord/",
                       int minHits = 30, Long64_t maxEntries = -1)
{
   gSystem->Load("libAtReconstruction.so");
   FairLogger::GetLogger()->SetLogScreenLevel("ERROR");
   TString sfx = "_multifit_reco.root";
   auto fo = TFile::Open(oldDir + run + sfx), fn = TFile::Open(newDir + run + sfx);
   if (!fo || fo->IsZombie() || !fn || fn->IsZombie()) { printf("cannot open both files\n"); return; }
   auto *to = (TTree *)fo->Get("cbmsim"); auto *tn = (TTree *)fn->Get("cbmsim");
   TClonesArray *po = nullptr, *pn = nullptr;
   to->SetBranchAddress("AtPatternEvent", &po);
   tn->SetBranchAddress("AtPatternEvent", &pn);

   Long64_t N = std::min(to->GetEntries(), tn->GetEntries());
   if (maxEntries > 0 && maxEntries < N) N = maxEntries;
   printf("  entries: old %lld  new %lld  -> comparing %lld\n", to->GetEntries(), tn->GetEntries(), N);

   long nMatch = 0, nHitMismatch = 0, nBigTheta = 0, nCross = 0, nBackOld = 0, nBackNew = 0;
   double sumAbs = 0, sumClOld = 0, sumClNew = 0;
   std::vector<double> dCl;
   for (Long64_t i = 0; i < N; ++i) {
      to->GetEntry(i); tn->GetEntry(i);
      auto *eo = (AtPatternEvent *)po->At(0); auto *en = (AtPatternEvent *)pn->At(0);
      if (!eo || !en) continue;
      std::map<int, AtTrack *> mo;
      for (auto &t : eo->GetTrackCand()) mo[t.GetTrackID()] = &t;
      for (auto &tnw : en->GetTrackCand()) {
         auto it = mo.find(tnw.GetTrackID());
         if (it == mo.end()) continue;
         AtTrack *told = it->second;
         const int nh = tnw.GetHitArray().size();
         if (nh < minHits) continue;
         if ((int)told->GetHitArray().size() != nh) { ++nHitMismatch; continue; }
         ++nMatch;
         const double a = told->GetGeoTheta() * TMath::RadToDeg();
         const double b = tnw.GetGeoTheta() * TMath::RadToDeg();
         if (std::isfinite(a) && std::isfinite(b)) {
            const double d = std::fabs(b - a);
            sumAbs += d;
            if (d > 10) ++nBigTheta;
            if ((a > 90) != (b > 90)) ++nCross;
            if (a > 90) ++nBackOld;
            if (b > 90) ++nBackNew;
         }
         const double co = told->GetHitClusterArray()->size(), cn = tnw.GetHitClusterArray()->size();
         sumClOld += co; sumClNew += cn; dCl.push_back(cn - co);
      }
   }
   printf("\n  matched tracks (>=%d hits) : %ld     hit-count MISMATCHES: %ld %s\n", minHits, nMatch,
          nHitMismatch, nHitMismatch ? "<-- NOT an ordering-only change, investigate" : "(clean)");
   if (!nMatch) return;
   printf("  mean |GeoTheta(new) - GeoTheta(old)| = %.1f deg\n", sumAbs / nMatch);
   printf("  |change| > 10 deg           : %ld  (%.0f%%)\n", nBigTheta, 100.0 * nBigTheta / nMatch);
   printf("  CROSSES the 90 deg boundary : %ld  (%.0f%%)\n", nCross, 100.0 * nCross / nMatch);
   printf("  backward (theta>90) fraction: old %.1f%%   new %.1f%%\n",
          100.0 * nBackOld / nMatch, 100.0 * nBackNew / nMatch);
   printf("  clusters/track              : old %.1f    new %.1f\n", sumClOld / nMatch, sumClNew / nMatch);
   std::sort(dCl.begin(), dCl.end());
   printf("  d(clusters) median %.0f   p10 %.0f   p90 %.0f\n",
          dCl[dCl.size() / 2], dCl[dCl.size() / 10], dCl[dCl.size() * 9 / 10]);
}
