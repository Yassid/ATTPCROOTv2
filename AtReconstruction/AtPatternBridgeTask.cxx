#include "AtPatternBridgeTask.h"

#include "AtPatternCircle2D.h"
#include "AtPatternEvent.h"
#include "AtTrack.h"
#include "AtTrackTransformer.h"

#include <FairLogger.h>
#include <FairRootManager.h>

#include <TClonesArray.h>

#include <cmath>

ClassImp(AtPatternBridgeTask);

InitStatus AtPatternBridgeTask::Init()
{
   auto *ioMan = FairRootManager::Instance();
   if (ioMan == nullptr) {
      LOG(error) << "AtPatternBridgeTask::Init: no FairRootManager";
      return kERROR;
   }
   fPatternEventArray = dynamic_cast<TClonesArray *>(ioMan->GetObject(fInputBranch.Data()));
   if (fPatternEventArray == nullptr) {
      LOG(error) << "AtPatternBridgeTask: branch '" << fInputBranch << "' not found";
      return kERROR;
   }
   LOG(info) << "AtPatternBridgeTask: reading from " << fInputBranch
             << "  Smooth3D r=" << fClusterRadius << " mm, d=" << fClusterDistance << " mm";
   return kSUCCESS;
}

void AtPatternBridgeTask::Exec(Option_t * /*option*/)
{
   if (fPatternEventArray->GetEntriesFast() == 0)
      return;
   auto *patternEvent = dynamic_cast<AtPatternEvent *>(fPatternEventArray->At(0));
   if (patternEvent == nullptr)
      return;

   AtTools::AtTrackTransformer transformer;
   auto &tracks = patternEvent->GetTrackCand();
   for (auto &tr : tracks) {
      // 1) Pull circle parameters out of the AtPatternCircle2D, if available.
      const auto *pat = tr.GetPattern();
      const auto *circ = dynamic_cast<const AtPatterns::AtPatternCircle2D *>(pat);
      if (circ != nullptr) {
         auto center = circ->GetCenter();
         double R = circ->GetRadius();
         tr.SetGeoCenter({center.X(), center.Y()});
         tr.SetGeoRadius(R);
      } else {
         continue; // Pattern not a circle — UKF will skip.
      }

      // 2) Cluster the hits along the track via Smooth3D.
      tr.ResetHitClusterArray();
      transformer.ClusterizeSmooth3D(tr, fClusterRadius, fClusterDistance);

      // 3) Estimate GeoTheta from the cluster Z-progression: the helix dip
      // angle. Use first-vs-last cluster 3D vector.
      auto *clusters = tr.GetHitClusterArray();
      if (clusters == nullptr || clusters->size() < 2) {
         tr.SetGeoTheta(M_PI / 2.); // fallback: 90 deg (track in pad plane)
         continue;
      }
      const auto &p0 = clusters->front().GetPosition();
      const auto &pN = clusters->back().GetPosition();
      const double dx = pN.X() - p0.X();
      const double dy = pN.Y() - p0.Y();
      const double dz = pN.Z() - p0.Z();
      const double r3 = std::sqrt(dx * dx + dy * dy + dz * dz);
      if (r3 < 1e-3) {
         tr.SetGeoTheta(M_PI / 2.);
         continue;
      }
      tr.SetGeoTheta(std::acos(std::abs(dz) / r3));
      tr.SetGeoPhi(std::atan2(dy, dx));
   }
}
