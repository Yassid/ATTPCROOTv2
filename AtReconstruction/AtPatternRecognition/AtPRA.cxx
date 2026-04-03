#include "AtPRA.h"

#include "AtContainerManip.h"
#include "AtHit.h"     // for AtHit, XYZPoint
#include "AtPattern.h" // for AtPattern
#include "AtPatternCircle2D.h"
#include "AtPatternEvent.h"
#include "AtPatternLine.h"
#include "AtPatternTypes.h" // for PatternType, PatternTy...
#include "AtSampleConsensus.h"
#include "AtTrack.h" // for XYZPoint, AtTrack

#include <FairLogger.h>

#include <Math/Point3D.h>     // for PositionVector3D, Cart...
#include <Math/Vector2D.h>    // for PositionVector3D, Cart...
#include <Math/Vector2Dfwd.h> // for XYVector
#include <Math/Vector3D.h>    // for DisplacementVector3D
#include <TGraph.h>           // for TGraph
#include <TMath.h>            // for Power, Sqrt, ATan2, Pi

#include <algorithm> // for max, for_each, copy_if
#include <cmath>     // for fabs, acos
#include <cstddef>   // for size_t
#include <exception> // for exception
#include <iostream>  // for operator<<, basic_ostream
#include <memory>    // for shared_ptr, __shared_p...
ClassImp(AtPATTERN::AtPRA);

/**
 * @brief Set initial parameters for HC.
 *
 * In track, sets GeoTheta, GeoPhi, GeoCenter, GeoRadius.
 */
void AtPATTERN::AtPRA::SetTrackInitialParameters(AtTrack &track)
{

   /*
std::cout << " Processing track with " << track.GetHitArray().size() << " points."
             << "\n";
   for (auto &hit : track.GetHitArray())
      std::cout << hit.GetTimeStamp() << " ";
   std::cout << std::endl;
   */

   // Select hits for the circle fit from the vertex end.
   // Use fRadiusFitFraction of hits from the vertex end (highest Z_digi)
   // where momentum is highest and the track is most circular.
   // Capped by fMinHitsRadius (floor) and fMaxHitsRadius (ceiling).
   auto &allHits = track.GetHitArray();
   int nTotal = allHits.size();
   int nFromFraction = std::max(3, static_cast<int>(nTotal * fRadiusFitFraction));
   int nHitsForFit = std::max(std::min(fMinHitsRadius, nTotal),
                              std::min(nFromFraction, std::min(fMaxHitsRadius, nTotal)));

   // Build the subset from the vertex end (highest Z_digi = last in array)
   std::vector<const AtHit *> hitsForFit;
   int startIdx = std::max(0, nTotal - nHitsForFit);
   for (int i = startIdx; i < nTotal; i++)
      hitsForFit.push_back(allHits[i].get());

   LOG(debug) << "AtPRA: circle fit using " << hitsForFit.size() << "/" << nTotal << " hits";

   SampleConsensus::AtSampleConsensus RansacSmoothRadius;
   RansacSmoothRadius.SetPatternType(AtPatterns::PatternType::kCircle2D);
   RansacSmoothRadius.SetMinHitsPattern(0.1 * hitsForFit.size());
   RansacSmoothRadius.SetDistanceThreshold(6.0);
   RansacSmoothRadius.SetNumIterations(1000);
   auto circularTracks = RansacSmoothRadius.Solve(hitsForFit).GetTrackCand();

   if (!circularTracks.empty()) {

      auto &hits = circularTracks.at(0).GetHitArray();

      // Fit circle on the selected subset
      auto circle = std::make_unique<AtPatterns::AtPatternCircle2D>();
      circle->AtPattern::FitPattern(hitsForFit);

      auto center = circle->GetCenter();
      auto radius = circle->GetRadius();

      track.SetGeoCenter({center.X(), center.Y()});
      track.SetGeoRadius(radius);
      track.SetPattern(circle->Clone());

      circularTracks.at(0).SetPattern(std::move(circle));

      // std::vector<double> wpca;
      std::vector<double> whit;
      std::vector<double> arclength;

      auto arclengthGraph = std::make_unique<TGraph>();

      auto posPCA = hits.at(0)->GetPosition();
      auto refPosOnCircle = posPCA - center;
      auto refAng = refPosOnCircle.Phi(); // Bounded between (-Pi,Pi]

      std::vector<AtHit> thetaHits;

      // The number of times we have crossed the -Y axis.
      // Increase by 1 when moving from +y to -y
      // Decrease by 1 when moving from -y to +y
      int numYCross = 0;
      int lastYSign = GetSign(refPosOnCircle.Y());

      for (size_t i = 0; i < hits.size(); ++i) {

         auto pos = hits.at(i)->GetPosition();
         auto posOnCircle = pos - center;
         auto angleHit = posOnCircle.Phi();

         // Check for a move over -Y axis
         int currYSign = GetSign(posOnCircle.Y());

         // If we have moved over the Y axis in some direction
         // last = (0 or -1) and current = 1 -> numCross--
         // last = (0 or 1) and current = -1 -> numCross++
         if (posOnCircle.X() < 0 && lastYSign != currYSign) {
            numYCross -= currYSign;
         }
         lastYSign = currYSign;
         angleHit += 2 * M_PI * numYCross;

         whit.push_back(angleHit);
         arclength.push_back((radius * (refAng - whit.at(i))));

         arclengthGraph->SetPoint(arclengthGraph->GetN(), arclength.at(i), pos.Z());

         if (track.GetTrackID() > -1)
            LOG(debug2) << posOnCircle.X() << "  " << posOnCircle.Y() << " " << pos.Z() << " "
                        << hits.at(i)->GetTimeStamp() << " " << arclength.back() << "\n";

         // Add a hit in the Arc legnth - Z plane
         Double_t xPos = arclength.at(i);
         Double_t yPos = pos.Z();
         Double_t zPos = i * 1E-19;

         thetaHits.emplace_back(i, hits.at(i)->GetPadNum(), XYZPoint(xPos, yPos, zPos), hits.at(i)->GetCharge());
      }

      // TF1 *f1 = new TF1("f1", "pol1", -500, 500);
      // TF1 * f1 = new TF1("f1",[](double *x, double *p) { return (p[0]+p[1]*x[0]); },-500,500,2);
      // TF1 * f1 = new TF1("f1","[0]+[1]*x",-500,500);
      // TF1 *f1 = new TF1("f1", fitf, -500, 500, 2);
      // arclengthGraph->Fit(f1, "R");
      // auto slope = ROOT::Math::XYVector(1, f1->GetParameter(1)).Unit();
      // std::cout << " Slope " << slope << "\n";

      Double_t angle = 0.0;
      Double_t phi0 = 0.0;

      try {

         if (thetaHits.size() > 0) {
            // std::cout<<" RANSAC Theta "<<"\n";
            std::vector<AtTrack> thetaTracks;

            SampleConsensus::AtSampleConsensus RansacTheta;
            RansacTheta.SetPatternType(AtPatterns::PatternType::kLine);
            RansacTheta.SetMinHitsPattern(0.1 * thetaHits.size());
            RansacTheta.SetDistanceThreshold(6.0);
            RansacTheta.SetFitPattern(true);
            thetaTracks = RansacTheta.Solve(thetaHits).GetTrackCand();

            if (thetaTracks.size() > 0) {

               // NB: Only the most intense line is taken, if any
               auto line = dynamic_cast<const AtPatterns::AtPatternLine *>(thetaTracks.at(0).GetPattern());
               LOG(info) << "Track ID: " << track.GetTrackID() << " with " << track.GetHitArray().size()
                         << " hits. Fit " << thetaTracks[0].GetHitArray().size() << "/" << thetaHits.size()
                         << " hits in first of " << thetaTracks.size() << " tracks";

               auto dirTheta = line->GetDirection();
               auto dir2D = ROOT::Math::XYVector(dirTheta.X(), dirTheta.Y()).Unit();

               int sign = 0;

               if (dir2D.X() * dir2D.Y() < 0)
                  sign = -1;
               else
                  sign = 1;

               if (dir2D.X() != 0) {
                  angle = acos(sign * fabs(dir2D.Y())) * TMath::RadToDeg();
               }

               LOG(info) << "Setting theta geo to: " << angle << " with ransac direction: " << dir2D;
               for (int i = 1; i < thetaTracks.size(); ++i) {
                  auto ranLine = dynamic_cast<const AtPatterns::AtPatternLine *>(thetaTracks.at(i).GetPattern());
                  auto ranDir = ROOT::Math::XYVector(ranLine->GetDirection().X(), ranLine->GetDirection().Y()).Unit();
                  LOG(info) << "RANSAC direction " << i << " with " << thetaTracks[i].GetHitArray().size() << "/"
                            << thetaHits.size() << ": " << ranDir;
               }

               auto temp2 = posPCA - center;
               phi0 = TMath::ATan2(temp2.Y(), temp2.X());

            } // thetaTracks
            track.SetGeoTheta(angle * TMath::Pi() / 180.0);
            track.SetGeoPhi(phi0);

         } // if
         else {
            LOG(info) << "Track ID: " << track.GetTrackID() << " with " << track.GetHitArray().size()
                      << " hits does not have enough theta hits to get angle";
         }

      } catch (std::exception &e) {

         std::cout << " AtPRA::SetTrackInitialParameters - Exception caught : " << e.what() << "\n";
      }

   } // end if (!circularTracks->empty())
}

Double_t fitf(Double_t *x, Double_t *par)
{

   return par[0] + par[1] * x[0];
}

void AtPATTERN::AtPRA::OrderClustersAlongTrack(AtTrack &track)
{
   auto *clusters = track.GetHitClusterArray();
   int nCl = clusters->size();
   if (nCl < 3)
      return;

   // Find seed: cluster with highest Z (vertex end, closest to beam entrance)
   int seedIdx = 0;
   for (int i = 1; i < nCl; i++) {
      if (clusters->at(i).GetPosition().Z() > clusters->at(seedIdx).GetPosition().Z())
         seedIdx = i;
   }

   // Nearest-neighbor walk from seed
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

   // Rebuild cluster array in the new order
   std::vector<AtHitCluster> ordered;
   ordered.reserve(order.size());
   for (int idx : order)
      ordered.push_back(clusters->at(idx));

   *clusters = std::move(ordered);
}

void AtPATTERN::AtPRA::MergeTrackFragments(std::vector<AtTrack> &tracks, double maxDist)
{
   if (tracks.size() <= 1)
      return;

   // For each pair of tracks, check if the endpoint of one is close to
   // the endpoint of another. If so, merge the smaller into the larger.
   bool merged = true;
   while (merged) {
      merged = false;
      for (size_t i = 0; i < tracks.size() && !merged; i++) {
         auto *clI = tracks[i].GetHitClusterArray();
         if (clI->empty())
            continue;
         auto posI_front = clI->front().GetPosition();
         auto posI_back = clI->back().GetPosition();

         for (size_t j = i + 1; j < tracks.size() && !merged; j++) {
            auto *clJ = tracks[j].GetHitClusterArray();
            if (clJ->empty())
               continue;
            auto posJ_front = clJ->front().GetPosition();
            auto posJ_back = clJ->back().GetPosition();

            // Check all 4 endpoint combinations
            double d1 = (posI_back - posJ_front).R();
            double d2 = (posI_front - posJ_back).R();
            double d3 = (posI_back - posJ_back).R();
            double d4 = (posI_front - posJ_front).R();
            double dMin = std::min({d1, d2, d3, d4});

            // Check circle consistency: fragments on the same physical track
            // should have similar center, radius, and theta from RANSAC.
            double radiusI = tracks[i].GetGeoRadius();
            double radiusJ = tracks[j].GetGeoRadius();
            auto centerI = tracks[i].GetGeoCenter();
            auto centerJ = tracks[j].GetGeoCenter();
            double thetaI = tracks[i].GetGeoTheta() * 180.0 / TMath::Pi();
            double thetaJ = tracks[j].GetGeoTheta() * 180.0 / TMath::Pi();

            // Relative radius difference
            double radiusDiff = (radiusI > 0 && radiusJ > 0)
                                   ? std::abs(radiusI - radiusJ) / std::max(radiusI, radiusJ) * 100.0
                                   : 999;
            // Center distance
            double centerDist = std::sqrt((centerI.first - centerJ.first) * (centerI.first - centerJ.first) +
                                          (centerI.second - centerJ.second) * (centerI.second - centerJ.second));
            // Theta difference
            double thetaDiff = std::abs(thetaI - thetaJ);

            // Merge if: endpoints close AND (similar circle OR similar theta)
            bool circleMatch = (radiusDiff < 30.0 && centerDist < 20.0); // within 30% radius, 20mm center
            bool thetaMatch = (thetaDiff < fMergeAngleThreshold);
            bool shouldMerge = (dMin < maxDist) && (circleMatch || thetaMatch);

            if (shouldMerge) {
               // Merge j into i: move all hits from j to i
               for (auto &hit : tracks[j].GetHitArray())
                  tracks[i].AddHit(hit->Clone());

               // Remove track j
               tracks.erase(tracks.begin() + j);

               // Re-cluster and re-order the merged track
               tracks[i].ResetHitClusterArray();
               fTrackTransformer->ClusterizeSmooth3D(tracks[i], fClusterRadius > 0 ? fClusterRadius : 10.0,
                                                     fClusterDistance > 0 ? fClusterDistance : 20.0);
               OrderClustersAlongTrack(tracks[i]);

               merged = true;
               LOG(info) << "Merged fragments: " << tracks[i].GetHitArray().size() << " hits"
                         << " (d=" << dMin << "mm, dR=" << radiusDiff << "%, dCenter=" << centerDist
                         << "mm, dTheta=" << thetaDiff << "deg)";
            }
         }
      }
   }
}

void AtPATTERN::AtPRA::PruneTrack(AtTrack &track)
{
   auto &hitArray = track.GetHitArray();

   std::cout << "    === Prunning track : " << track.GetTrackID() << "\n";
   std::cout << "      = Hit Array size : " << hitArray.size() << "\n";

   for (auto iHit = 0; iHit < hitArray.size(); ++iHit) {

      try {
         bool isNoise = kNN(hitArray, *hitArray.at(iHit), fKNN); // Returns true if hit is an outlier

         if (isNoise) {
            // std::cout<<" Hit "<<iHit<<" flagged as outlier. "<<"\n";
            hitArray.erase(hitArray.begin() + iHit);
         }
      } catch (std::exception &e) {

         std::cout << " AtPRA::PruneTrack - Exception caught : " << e.what() << "\n";
      }
   }

   std::cout << "      = Hit Array size after prunning : " << hitArray.size() << "\n";
}

bool AtPATTERN::AtPRA::kNN(const std::vector<std::unique_ptr<AtHit>> &hits, AtHit &hitRef, int k)
{

   std::vector<Double_t> distances;
   distances.reserve(hits.size());

   std::for_each(hits.begin(), hits.end(), [&distances, &hitRef](const std::unique_ptr<AtHit> &hit) {
      distances.push_back(TMath::Sqrt((hitRef.GetPosition() - hit->GetPosition()).Mag2()));
   });

   std::sort(distances.begin(), distances.end(), [](Double_t a, Double_t b) { return a < b; });

   Double_t mean = 0.0;
   Double_t stdDev = 0.0;

   if (k > hits.size())
      k = hits.size();

   // Compute mean distance of kNN
   for (auto i = 0; i < k; ++i)
      mean += distances.at(i);

   mean /= k;

   // Compute std dev
   for (auto i = 0; i < k; ++i)
      stdDev += TMath::Power((distances.at(i) - mean), 2);

   stdDev = TMath::Sqrt(stdDev / k);

   // Compute threshold
   Double_t T = mean + stdDev * fStdDevMulkNN;

   return (T < fkNNDist) ? false : true;
}

/*
void AtPATTERN::AtPRA::SetTrackCurvature(AtTrack &track)
{
   std::cout << " new track  "
             << "\n";
   std::vector<double> radius_vec;
   std::vector<AtHit> hitArray = track.GetHitArray();
   int nstep = 0.60 * hitArray.size(); // 20% of the hits to calculate the radius of curvature with less fluctuations

   for (Int_t iHit = 0; iHit < (hitArray.size() - nstep); iHit++) {

      AtHit hitA = hitArray.at(iHit);
      AtHit hitB = hitArray.at((int)(iHit + (nstep / 2.0)));
      AtHit hitC = hitArray.at((int)(iHit + nstep));

      // std::cout<<nstep<<" "<<iHit<<"  "<<(int)(iHit+(nstep/2.0))<<"  "<<(int)(iHit+nstep)<<"\n";

      auto posA = hitA.GetPosition();
      auto posB = hitB.GetPosition();
      auto posC = hitC.GetPosition();

      double slopeAB = (posB.Y() - posA.Y()) / (posB.X() - posA.X()); // m1
      double slopeBC = (posC.Y() - posB.Y()) / (posC.X() - posB.X()); // m2

      double centerX = (slopeAB * slopeBC * (posA.Y() - posC.Y()) + slopeBC * (posB.X() + posA.X()) -
                        slopeAB * (posB.X() + posC.X())) /
                       (2.0 * (slopeBC - slopeAB));

      double centerY = (-1 / slopeAB) * (centerX - (posB.X() + posA.X()) / 2.0) + (posB.Y() + posA.Y()) / 2.0;

      // std::cout<<" Center "<<centerX<<" - "<<centerY<<"\n";

      double radiusA = TMath::Sqrt(TMath::Power(posA.X() - centerX, 2) + TMath::Power(posA.Y() - centerY, 2));
      radius_vec.push_back(radiusA);
      double radiusB = TMath::Sqrt(TMath::Power(posB.X() - centerX, 2) + TMath::Power(posB.Y() - centerY, 2));
      radius_vec.push_back(radiusB);
      double radiusC = TMath::Sqrt(TMath::Power(posC.X() - centerX, 2) + TMath::Power(posC.Y() - centerY, 2));
      radius_vec.push_back(radiusC);
   }
}
*/
