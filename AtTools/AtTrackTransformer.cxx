#include "AtTrackTransformer.h"
// IWYU pragma: no_include <ext/alloc_traits.h>

#include "AtHit.h"        // for AtHit, AtHit::XYZPoint
#include "AtHitCluster.h" // for AtHitCluster
#include "AtTrack.h"      // for XYZPoint, AtTrack

#include <Math/Point3D.h> // for PositionVector3D, Cart...
#include <Math/Point3Dfwd.h>
#include <Math/Vector3D.h>  // for DisplacementVector3D
#include <TMath.h>          // for Power, Sqrt, ATan2, Pi
#include <TMatrixDSymfwd.h> // for TMatrixDSym
#include <TMatrixTSym.h>    // for TMatrixTSym
#include <TVector3.h>       // for TVector3

#include <algorithm> // for max, for_each, copy_if
#include <iterator>  // for back_insert_iterator
#include <memory>    // for shared_ptr, __shared_p...
#include <vector>    // for vector

AtTools::AtTrackTransformer::AtTrackTransformer() = default;
AtTools::AtTrackTransformer::~AtTrackTransformer() = default;

void AtTools::AtTrackTransformer::SetDiffusionParams(double coefT, double coefL, double driftVel, double tbTime,
                                                     double padResXY)
{
   fCoefT = coefT;
   fCoefL = coefL;
   fDriftVel = driftVel;
   fTBTime = tbTime;
   if (padResXY > 0)
      fPadResXY = padResXY;
}
using XYZPoint = ROOT::Math::XYZPoint;

void AtTools::AtTrackTransformer::ClusterizeSmooth3D(AtTrack &track, Float_t radius, Float_t distance)
{
   std::vector<AtHit> hitArray = track.GetHitArrayObject();
   std::vector<AtHit> hitTBArray;
   int clusterID = 0;

   // std::cout<<" ================================================================= "<<"\n";
   // std::cout<<" Clusterizing track : "<<track.GetTrackID()<<"\n";

   /*for(auto iHits=0;iHits<hitArray->size();++iHits)
     {
       TVector3 pos    = hitArray->at(iHits).GetPosition();
       double Q = hitArray->at(iHits).GetCharge();
       int TB          = hitArray->at(iHits).GetTimeStamp();
       //std::cout<<" Pos : "<<pos.X()<<" - "<<pos.Y()<<" - "<<pos.Z()<<" - TB : "<<TB<<" - Charge : "<<Q<<"\n";
       }*/

   // Diffusion and resolution parameters.
   // The transverse/longitudinal diffusion sigmas (in mm) at drift distance z_mm are:
   //   σ_T(z) = 10 * sqrt(CoefT * 2 * z_mm / (10 * vDrift))   [mm]
   //   σ_L(z) = 10 * sqrt(CoefL * 2 * z_mm / (10 * vDrift))   [mm]
   // where CoefT/L are in cm²/µs, vDrift in cm/µs, and z_mm is the drift distance in mm.
   // This matches the digitization in AtClusterize::getTransverseDiffusion/getLongitudinalDiffusion.
   Double_t driftVel = fDriftVel;       // cm/us
   Double_t samplingRate = fTBTime;     // us
   Double_t padResXY = fPadResXY;       // mm (transverse pad resolution)
   Double_t padResZ = fPadResXY * 1.5;  // mm (longitudinal — pads are typically longer in Z)

   if (hitArray.size() > 0) {

      auto refPos = hitArray.at(0).GetPosition(); // First hit
      // TODO: Create a clustered hit from the very first hit (test)

      for (auto iHit = 0; iHit < hitArray.size(); ++iHit) {

         auto hit = hitArray.at(iHit);

         // Check distance with respect to reference Hit
         Double_t distRef = TMath::Sqrt((hit.GetPosition() - refPos).Mag2());

         if (distRef < distance) {

            continue;

         } else {

            // std::cout<<" Clustering "<<iHit<<" of "<<hitArray->size()<<"\n";
            // std::cout<<" Distance to reference : "<<distRef<<"\n";
            // std::cout<<" Reference position : "<<refPos.X()<<" - "<<refPos.Y()<<" - "<<refPos.Z()<<" -
            // "<<refPos.Mag()<<"\n";

            Double_t clusterQ = 0.0;
            hitTBArray.clear();
            std::copy_if(
               hitArray.begin(), hitArray.end(), std::back_inserter(hitTBArray),
               [&refPos, radius](AtHit &hitIn) { return TMath::Sqrt((hitIn.GetPosition() - refPos).Mag2()) < radius; });

            // std::cout<<" Clustered "<<hitTBArray.size()<<" Hits "<<"\n";

            if (hitTBArray.size() > 0) {
               double x = 0, y = 0, z = 0;
               double var_x = 0, var_y = 0, var_z = 0; // accumulate variance (σ²)

               int timeStamp = 0;
               std::shared_ptr<AtHitCluster> hitCluster = std::make_shared<AtHitCluster>();
               hitCluster->SetClusterID(clusterID);
               Double_t hitQ = 0.0;
               Double_t hitQ2 = 0.0; // sum of charge² for weighted variance
               int nHits = 0;

               for (auto &hitInQ : hitTBArray) {
                  XYZPoint pos = hitInQ.GetPosition();
                  double q = hitInQ.GetCharge();
                  x += pos.X() * q;
                  y += pos.Y() * q;
                  z += pos.Z();
                  hitQ += q;
                  hitQ2 += q * q;
                  timeStamp += hitInQ.GetTimeStamp();
                  nHits++;

                  // Per-hit position variance (σ² in mm²).
                  // Drift distance: pos.Z() is in mm (digi frame, distance from pad plane).
                  // Convert to drift time: t_drift = z_mm / (10 * vDrift_cm_per_us) [µs]
                  double driftTime = pos.Z() / (10.0 * driftVel); // µs

                  // Transverse diffusion variance: σ_T² = (10*sqrt(CoefT*2*t))² = 100*CoefT*2*t [mm²]
                  double varT = 100.0 * fCoefT * 2.0 * driftTime; // mm²

                  // Longitudinal diffusion: σ_L in time [µs] → convert to mm
                  // σ_L_us = sqrt(CoefL*2*t), σ_L_mm = σ_L_us * vDrift * 10
                  double varL = 100.0 * fCoefL * 2.0 * driftTime; // mm²

                  // Time bucket resolution: uniform distribution → σ² = (vDrift*tbTime*10)²/12 [mm²]
                  double tbRes_mm = driftVel * samplingRate * 10.0; // mm per time bucket
                  double varTB = tbRes_mm * tbRes_mm / 12.0;       // mm²

                  // Charge-weighted accumulation of per-hit variance
                  var_x += q * q * (padResXY * padResXY + varT);
                  var_y += q * q * (padResXY * padResXY + varT);
                  var_z += q * q * (padResZ * padResZ + varTB + varL);
               }

               // Charge-weighted mean position
               x /= hitQ;
               y /= hitQ;
               z /= nHits;
               timeStamp /= nHits;

               // Variance of the charge-weighted mean: Var(x_bar) = Σ(w²*σ²) / (Σw)²
               // where w = charge. This gives the 1/√N improvement for the cluster.
               double sigma_x = TMath::Sqrt(var_x) / hitQ;
               double sigma_y = TMath::Sqrt(var_y) / hitQ;
               double sigma_z = TMath::Sqrt(var_z) / hitQ;

               XYZPoint clustPos(x, y, z);
               Bool_t checkDistance = kTRUE;

               // Check distance with respect to existing clusters
               for (auto iClusterHit : *track.GetHitClusterArray()) {
                  if (TMath::Sqrt((iClusterHit.GetPosition() - clustPos).Mag2()) < distance) {
                     // std::cout<<" Cluster with less than  : "<<distance<<" found "<<"\n";
                     checkDistance = kFALSE;
                     continue;
                  }
               }

               if (checkDistance) {
                  hitCluster->SetCharge(hitQ);
                  hitCluster->SetPosition({x, y, z});
                  hitCluster->SetTimeStamp(timeStamp);
                  TMatrixDSym cov(3); // TODO: Setting covariant matrix based on pad size and drift time resolution.
                                      // Using estimations for the moment.
                  cov(0, 1) = 0;
                  cov(1, 2) = 0;
                  cov(2, 0) = 0;
                  cov(0, 0) = TMath::Power(sigma_x, 2); // 0.04;
                  cov(1, 1) = TMath::Power(sigma_y, 2); // 0.04;
                  cov(2, 2) = TMath::Power(sigma_z, 2); // 0.01;
                  hitCluster->SetCovMatrix(cov);
                  ++clusterID;
                  track.AddClusterHit(hitCluster);
               }
            }
         }

         // Sanity check
         /*std::cout<<" Hits for cluster "<<iHit<<" centered in "<<refPos.X()<<" - "<<refPos.Y()<<"-"<<refPos.Z()<<"\n";
    for(auto iHits=0;iHits<hitTBArray.size();++iHits)
         {
           TVector3 pos    = hitTBArray.at(iHits).GetPosition();
           double Q = hitTBArray.at(iHits).GetCharge();
           int TB          = hitTBArray.at(iHits).GetTimeStamp();
           std::cout<<" Pos : "<<pos.X()<<" - "<<pos.Y()<<" - "<<pos.Z()<<" - TB : "<<TB<<" - Charge : "<<Q<<"\n";
      std::cout<<" Distance to cluster center "<<TMath::Abs((track.GetHitClusterArray()->back().GetPosition() -
    pos).Mag())<<"\n";
    }
         std::cout<<"=================================================="<<"\n";*/

         refPos = hitArray.at(iHit).GetPosition();

         //} // if distance

      } // for

      // Smoothing track
      std::vector<AtHitCluster> *hitClusterArray = track.GetHitClusterArray();
      radius /= 2.0;
      std::vector<std::shared_ptr<AtHitCluster>> hitClusterBuffer;

      // std::cout<<" Hit cluster array size "<<hitClusterArray->size()<<"\n";

      if (hitClusterArray->size() > 2) {

         for (auto iHitCluster = 0; iHitCluster < hitClusterArray->size() - 1;
              ++iHitCluster) // Calculating distances between pairs of clusters
         {

            XYZPoint clusBack = hitClusterArray->at(iHitCluster).GetPosition();
            XYZPoint clusForw = hitClusterArray->at(iHitCluster + 1).GetPosition();
            XYZPoint clusMidPos = clusBack + (clusForw - clusBack) * 0.5;
            std::vector<XYZPoint> renormClus{clusBack, clusMidPos};

            if (iHitCluster == (hitClusterArray->size() - 2))
               renormClus.push_back(clusForw);

            // Create a new cluster and renormalize the charge of the other with half the radius.
            for (auto iClus : renormClus) {
               hitTBArray.clear();
               std::copy_if(hitArray.begin(), hitArray.end(), std::back_inserter(hitTBArray),
                            [&iClus, radius](AtHit &hitIn) {
                               return TMath::Sqrt((hitIn.GetPosition() - iClus).Mag2()) < radius;
                            });

               if (hitTBArray.size() > 0) {
                  double x = 0, y = 0, z = 0;
                  double var_x = 0, var_y = 0, var_z = 0;

                  int timeStamp = 0;
                  std::shared_ptr<AtHitCluster> hitCluster = std::make_shared<AtHitCluster>();
                  hitCluster->SetClusterID(clusterID);
                  Double_t hitQ = 0.0;
                  int nHits = 0;

                  for (auto &hitInQ : hitTBArray) {
                     auto pos = hitInQ.GetPosition();
                     double q = hitInQ.GetCharge();
                     x += pos.X() * q;
                     y += pos.Y() * q;
                     z += pos.Z();
                     hitQ += q;
                     timeStamp += hitInQ.GetTimeStamp();
                     nHits++;

                     double driftTime = pos.Z() / (10.0 * driftVel);
                     double varT = 100.0 * fCoefT * 2.0 * driftTime;
                     double varL = 100.0 * fCoefL * 2.0 * driftTime;
                     double tbRes_mm = driftVel * samplingRate * 10.0;
                     double varTB = tbRes_mm * tbRes_mm / 12.0;

                     var_x += q * q * (padResXY * padResXY + varT);
                     var_y += q * q * (padResXY * padResXY + varT);
                     var_z += q * q * (padResZ * padResZ + varTB + varL);
                  }

                  x /= hitQ;
                  y /= hitQ;
                  z /= nHits;
                  timeStamp /= nHits;

                  double sigma_x = TMath::Sqrt(var_x) / hitQ;
                  double sigma_y = TMath::Sqrt(var_y) / hitQ;
                  double sigma_z = TMath::Sqrt(var_z) / hitQ;

                  TVector3 clustPos(x, y, z);
                  hitCluster->SetCharge(hitQ);
                  hitCluster->SetPosition({x, y, z});
                  hitCluster->SetTimeStamp(timeStamp);
                  TMatrixDSym cov(3); // TODO: Setting covariant matrix based on pad size and drift time resolution.
                                      // Using estimations for the moment.
                  cov(0, 1) = 0;
                  cov(1, 2) = 0;
                  cov(2, 0) = 0;
                  cov(0, 0) = TMath::Power(sigma_x, 2); // 0.04;
                  cov(1, 1) = TMath::Power(sigma_y, 2); // 0.04;
                  cov(2, 2) = TMath::Power(sigma_z, 2); // 0.01;
                  hitCluster->SetCovMatrix(cov);
                  ++clusterID;
                  hitClusterBuffer.push_back(hitCluster);

               } // hitTBArray size

            } // for iClus

         } // for HitArray

         // Remove previous clusters
         track.ResetHitClusterArray();

         // Adding new clusters
         for (auto iHitClusterRe : hitClusterBuffer) {

            track.AddClusterHit(iHitClusterRe);
         }

      } // Cluster array size

   } // if array size
}

void AtTools::AtTrackTransformer::ClusterizeByGroup(AtTrack &track, int hitsPerCluster)
{
   auto hitArray = track.GetHitArrayObject(); // copy as value (vector<AtHit>)
   if (hitArray.empty())
      return;

   // Diffusion/resolution params for covariance
   Double_t driftVel = fDriftVel;
   Double_t samplingRate = fTBTime;
   Double_t padResXY = fPadResXY;
   Double_t padResZ = fPadResXY * 1.5;

   int nHits = hitArray.size();
   int clusterID = 0;

   for (int start = 0; start < nHits; start += hitsPerCluster) {
      int end = std::min(start + hitsPerCluster, nHits);
      int groupSize = end - start;
      if (groupSize < 2)
         continue;

      // Charge-weighted centroid
      double x = 0, y = 0, z = 0;
      double totalQ = 0;
      double var_x = 0, var_y = 0, var_z = 0;
      int timeStamp = 0;

      for (int i = start; i < end; i++) {
         auto &hit = hitArray[i];
         auto pos = hit.GetPosition();
         double q = hit.GetCharge();

         x += pos.X() * q;
         y += pos.Y() * q;
         z += pos.Z();
         totalQ += q;
         timeStamp += hit.GetTimeStamp();

         // Per-hit variance
         double driftTime = pos.Z() / (10.0 * driftVel);
         double varT = 100.0 * fCoefT * 2.0 * driftTime;
         double varL = 100.0 * fCoefL * 2.0 * driftTime;
         double tbRes_mm = driftVel * samplingRate * 10.0;
         double varTB = tbRes_mm * tbRes_mm / 12.0;

         var_x += q * q * (padResXY * padResXY + varT);
         var_y += q * q * (padResXY * padResXY + varT);
         var_z += q * q * (padResZ * padResZ + varTB + varL);
      }

      if (totalQ <= 0)
         continue;

      x /= totalQ;
      y /= totalQ;
      z /= groupSize;
      timeStamp /= groupSize;

      double sigma_x = TMath::Sqrt(var_x) / totalQ;
      double sigma_y = TMath::Sqrt(var_y) / totalQ;
      double sigma_z = TMath::Sqrt(var_z) / totalQ;

      auto hitCluster = std::make_shared<AtHitCluster>();
      hitCluster->SetClusterID(clusterID++);
      hitCluster->SetCharge(totalQ);
      hitCluster->SetPosition({x, y, z});
      hitCluster->SetTimeStamp(timeStamp);

      TMatrixDSym cov(3);
      cov(0, 0) = sigma_x * sigma_x;
      cov(1, 1) = sigma_y * sigma_y;
      cov(2, 2) = sigma_z * sigma_z;
      cov(0, 1) = cov(1, 0) = 0;
      cov(0, 2) = cov(2, 0) = 0;
      cov(1, 2) = cov(2, 1) = 0;
      hitCluster->SetCovMatrix(cov);

      track.AddClusterHit(hitCluster);
   }
}

Bool_t AtTools::AtTrackTransformer::FindVertexTrack(AtTrack *trA, AtTrack *trB)
{
   // Determination of first hit distance. NB: Assuming both tracks have the same angle sign
   Double_t vertexA = 0.0;
   Double_t vertexB = 0.0;
   if (trA->GetGeoTheta() * TMath::RadToDeg() < 90) {
      auto iniClusterA = trA->GetHitClusterArray()->back();
      auto iniClusterB = trB->GetHitClusterArray()->back();
      vertexA = 1000.0 - iniClusterA.GetPosition().Z();
      vertexB = 1000.0 - iniClusterB.GetPosition().Z();
   } else if (trA->GetGeoTheta() * TMath::RadToDeg() > 90) {
      auto iniClusterA = trA->GetHitClusterArray()->front();
      auto iniClusterB = trB->GetHitClusterArray()->front();
      vertexA = iniClusterA.GetPosition().Z();
      vertexB = iniClusterB.GetPosition().Z();
   }

   return vertexA < vertexB;
}

const std::tuple<Double_t, Double_t> AtTools::AtTrackTransformer::GetPIDFromHits(AtTrack &track, Double_t theta)
{

   Double_t dedx = 0.0;
   Double_t eloss = 0.0;

   auto hitArray = &track.GetHitArray();
   std::size_t cnt = 0;

   if (theta < 90) {

      auto it = hitArray->rbegin();
      while (it != hitArray->rend()) {

         if (((Float_t)cnt / (Float_t)hitArray->size()) > 0.8)
            break;

         eloss += (*it).get()->GetCharge();
         dedx += (*it).get()->GetCharge();
         // std::cout<<(*it).GetCharge()<<"\n";
         it++;
         ++cnt;
      }
   } else if (theta > 90) {

      eloss += hitArray->at(0).get()->GetCharge();

      cnt = 1;
      for (auto iHitClus = 1; iHitClus < hitArray->size(); ++iHitClus) {

         if (((Float_t)cnt / (Float_t)hitArray->size()) > 0.8)
            break;

         eloss += hitArray->at(iHitClus).get()->GetCharge();
         dedx += hitArray->at(iHitClus).get()->GetCharge();
         // std::cout<<len<<" - "<<eloss<<" - "<<hitClusterArray->at(iHitClus).GetCharge()<<"\n";
         ++cnt;
      }
   }

   eloss /= cnt;

   return std::forward_as_tuple(dedx, eloss);
}

Bool_t AtTools::AtTrackTransformer::MergeTracks(std::vector<AtTrack *> *trackCandSource,
                                                std::vector<AtTrack> *trackDest, Bool_t enableSingleVertexTrack,
                                                Double_t clusterRadius, Double_t clusterDistance)
{

   Bool_t toMerge = kFALSE;

   Int_t addHitCnt = 0;
   // Find the track closer to vertex
   std::sort(trackCandSource->begin(), trackCandSource->end(),
             [this](AtTrack *trA, AtTrack *trB) { return FindVertexTrack(trA, trB); });

   // Track stitching from vertex
   AtTrack *vertexTrack = *trackCandSource->begin();

   if (enableSingleVertexTrack) {

      // Mark all tracks as merged
      for (auto track : *trackCandSource)
         track->SetIsMerged(kTRUE);

      trackDest->push_back(*vertexTrack);
      return true;
   }

   // Check if the candidate vertex track was merged
   if (vertexTrack->GetIsMerged())
      return kFALSE;
   else
      vertexTrack->SetIsMerged(kTRUE);

   // If enabled, choose only the track closest to vertex (i.e. first one of the collection of candidates)
   // TODO: Select by number of points

   for (auto it = trackCandSource->begin() + 1; it != trackCandSource->end(); ++it) {
      // NB: These tracks were previously marked to merge. If merging fails they should be discarded.
      AtTrack *trackToMerge = *(it);
      toMerge = kFALSE;

      // Skip trackes flagged as merged
      if (!trackToMerge->GetIsMerged()) {
         trackToMerge->SetIsMerged(kTRUE);
      } else
         continue;

      Double_t endVertexZ = 0.0;
      Double_t iniMergeZ = 0.0;
      std::cout << " Trying to merge ... "
                << "\n";
      std::cout << " Vertex track " << vertexTrack->GetTrackID() << " - Track to Merge " << trackToMerge->GetTrackID()
                << "\n";
      // Check relative position between end and begin of each track using Hit Clusters
      std::cout << " Vertex angle " << vertexTrack->GetGeoTheta() * TMath::RadToDeg() << "\n";
      if (vertexTrack->GetGeoTheta() * TMath::RadToDeg() < 90) {
         auto endClusterVertex = vertexTrack->GetHitClusterArray()->front();
         auto iniClusterMerge = trackToMerge->GetHitClusterArray()->back();
         // Check separation and relative distance
         endVertexZ = 1000.0 - endClusterVertex.GetPosition().Z();
         iniMergeZ = 1000.0 - iniClusterMerge.GetPosition().Z();

         Double_t distance = std::sqrt((iniClusterMerge.GetPosition() - endClusterVertex.GetPosition()).Mag2());
         // std::cout << " Distance between tracks " << distance << "\n";
         // std::cout << " Ini Merge " << iniMergeZ << " - endVertexZ " << endVertexZ << "\n";
         if (((iniMergeZ + 10.0) > endVertexZ) && distance < 200) {
            toMerge = kTRUE;
         }

      } else if (vertexTrack->GetGeoTheta() * TMath::RadToDeg() > 90) {
         auto endClusterVertex = vertexTrack->GetHitClusterArray()->back();
         auto iniClusterMerge = trackToMerge->GetHitClusterArray()->front();
         // Check separation and relative distance
         endVertexZ = endClusterVertex.GetPosition().Z();
         iniMergeZ = iniClusterMerge.GetPosition().Z();
         Double_t distance = std::sqrt((iniClusterMerge.GetPosition() - endClusterVertex.GetPosition()).Mag2());
         // std::cout<<" Distance between tracks "<<distance<<"\n";
         // std::cout<<" Ini Merge "<<iniMergeZ<<" - endVertexZ "<<endVertexZ<<"\n";
         if (((iniMergeZ + 10.0) > endVertexZ) &&
             distance < 100) { // NB: Distance between parts of the backward tracks is more critical
            toMerge = kTRUE;
         }
      }

      if (toMerge) {

         std::cout << " --- Merging Succeeded! Vertex track " << vertexTrack->GetTrackID() << " - Track to Merge "
                   << trackToMerge->GetTrackID() << "\n";
         for (const auto &hit : trackToMerge->GetHitArray()) {

            vertexTrack->AddHit(hit->Clone()); // TODO: Look at code and see if this can be a move instead of a copy
            ++addHitCnt;
         }

         // Reclusterize after merging
         vertexTrack->SortHitArrayTime();
         vertexTrack->ResetHitClusterArray();
         ClusterizeSmooth3D(
            *vertexTrack, clusterRadius,
            clusterDistance); // NB: It can be removed if we force reclusterization for any track in the mina program

         // TODO: Check if phi recalculatio is needed

      } else {
         std::cout << " --- Merging Failed ! Vertex track " << vertexTrack->GetTrackID() << " - Track to Merge "
                   << trackToMerge->GetTrackID() << "\n";
      }
   }

   trackDest->push_back(*vertexTrack);

   return toMerge;
}
