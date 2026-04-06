#include "AtTrackRefiner.h"

#include "AtHit.h"
#include "AtHitCluster.h"
#include "AtTrack.h"
#include "AtTrackTransformer.h"

#include <FairLogger.h>

#include <Math/Point3D.h>

#include <algorithm>
#include <cmath>
#include <memory>
#include <vector>

using XYZPoint = ROOT::Math::XYZPoint;

void AtPATTERN::AtTrackRefiner::OrderClustersAlongTrack(AtTrack &track)
{
   auto *clusters = track.GetHitClusterArray();
   int nCl = clusters->size();
   if (nCl < 3)
      return;

   int seedIdx = 0;
   for (int i = 1; i < nCl; i++) {
      if (clusters->at(i).GetPosition().Z() > clusters->at(seedIdx).GetPosition().Z())
         seedIdx = i;
   }

   std::vector<int> order;
   std::vector<bool> used(nCl, false);
   order.push_back(seedIdx);
   used[seedIdx] = true;
   for (int step = 1; step < nCl; step++) {
      auto current = clusters->at(order.back()).GetPosition();
      double bestDist = 1e9;
      int bestIdx = -1;
      for (int j = 0; j < nCl; j++) {
         if (used[j])
            continue;
         double d = (clusters->at(j).GetPosition() - current).R();
         if (d < bestDist) {
            bestDist = d;
            bestIdx = j;
         }
      }
      if (bestIdx < 0)
         break;
      order.push_back(bestIdx);
      used[bestIdx] = true;
   }

   std::vector<AtHitCluster> ordered;
   ordered.reserve(order.size());
   for (int idx : order)
      ordered.push_back(clusters->at(idx));

   *clusters = std::move(ordered);
}

void AtPATTERN::AtTrackRefiner::SelectAndMergeTracks(std::vector<AtTrack> &tracks,
                                                     AtTools::AtTrackTransformer &trackTransformer,
                                                     double clusterRadius, double clusterDistance,
                                                     double vertexRadiusXY, double mergeDist, double minLabTheta)
{
   if (tracks.empty())
      return;

   auto getVertexEnd = [](AtTrack &tr) -> XYZPoint {
      auto *cl = tr.GetHitClusterArray();
      if (cl->empty())
         return XYZPoint(0, 0, 0);
      return cl->front().GetPosition();
   };

   auto getFarEnd = [](AtTrack &tr) -> XYZPoint {
      auto *cl = tr.GetHitClusterArray();
      if (cl->empty())
         return XYZPoint(0, 0, 0);
      return cl->back().GetPosition();
   };

   auto xyDist = [](const XYZPoint &p) { return std::sqrt(p.X() * p.X() + p.Y() * p.Y()); };

   tracks.erase(std::remove_if(tracks.begin(), tracks.end(),
                               [minLabTheta, &getVertexEnd, &getFarEnd](AtTrack &tr) {
                                  auto vtx = getVertexEnd(tr);
                                  auto far = getFarEnd(tr);
                                  auto dir = far - vtx;
                                  if (dir.R() < 1e-3)
                                     return true;
                                  double cosThDigi = -dir.Z() / dir.R();
                                  double thDigi = std::acos(std::min(1.0, std::max(-1.0, cosThDigi))) * 180.0 / M_PI;
                                  double thLab = 180.0 - thDigi;
                                  return (thLab < minLabTheta || thLab > (180.0 - minLabTheta));
                               }),
                tracks.end());

   if (tracks.empty())
      return;

   double maxZ = -1e9;
   for (auto &tr : tracks) {
      auto vtxEnd = getVertexEnd(tr);
      if (vtxEnd.Z() > maxZ)
         maxZ = vtxEnd.Z();
   }

   std::vector<bool> isPrimary(tracks.size(), false);
   for (size_t i = 0; i < tracks.size(); i++) {
      auto vtxEnd = getVertexEnd(tracks[i]);
      double rXY = xyDist(vtxEnd);
      double dZ = std::abs(vtxEnd.Z() - maxZ);
      if (rXY < vertexRadiusXY && dZ < 50.0)
         isPrimary[i] = true;
   }

   int nPrimary = std::count(isPrimary.begin(), isPrimary.end(), true);
   LOG(info) << "SelectAndMerge: " << tracks.size() << " tracks, " << nPrimary << " primary (vertexR<"
             << vertexRadiusXY << "mm)";

   bool mergedAny = true;
   while (mergedAny) {
      mergedAny = false;
      for (size_t i = 0; i < tracks.size() && !mergedAny; i++) {
         if (!isPrimary[i])
            continue;

         auto farEnd = getFarEnd(tracks[i]);

         for (size_t j = 0; j < tracks.size() && !mergedAny; j++) {
            if (i == j || isPrimary[j])
               continue;

            auto *clJ = tracks[j].GetHitClusterArray();
            if (clJ->empty())
               continue;

            auto posJFront = clJ->front().GetPosition();
            auto posJBack = clJ->back().GetPosition();
            double d1 = (farEnd - posJFront).R();
            double d2 = (farEnd - posJBack).R();
            double dMin = std::min(d1, d2);

            if (dMin < mergeDist) {
               for (auto &hit : tracks[j].GetHitArray())
                  tracks[i].AddHit(hit->Clone());

               tracks.erase(tracks.begin() + j);
               isPrimary.erase(isPrimary.begin() + j);

               tracks[i].ResetHitClusterArray();
               trackTransformer.ClusterizeSmooth3D(tracks[i], clusterRadius > 0 ? clusterRadius : 10.0,
                                                   clusterDistance > 0 ? clusterDistance : 20.0);
               OrderClustersAlongTrack(tracks[i]);

               mergedAny = true;
               LOG(info) << "Merged fragment into primary: " << tracks[i].GetHitArray().size() << " hits (d=" << dMin
                         << "mm)";
            }
         }
      }
   }

   for (int i = tracks.size() - 1; i >= 0; i--) {
      if (!isPrimary[i]) {
         LOG(debug) << "Rejected isolated track: " << tracks[i].GetHitArray().size() << " hits";
         tracks.erase(tracks.begin() + i);
         isPrimary.erase(isPrimary.begin() + i);
      }
   }

   LOG(info) << "SelectAndMerge: " << tracks.size() << " tracks after selection";
}
