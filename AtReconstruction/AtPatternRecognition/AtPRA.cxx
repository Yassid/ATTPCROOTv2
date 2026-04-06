#include "AtPRA.h"

#include "AtContainerManip.h"
#include "AtHit.h"     // for AtHit, XYZPoint
#include "AtPattern.h" // for AtPattern
#include "AtPatternEvent.h"
#include "AtTrack.h" // for XYZPoint, AtTrack

#include <FairLogger.h>

#include <Math/Point3D.h>  // for PositionVector3D, Cart...
#include <Math/Vector3D.h> // for DisplacementVector3D
#include <TMath.h>         // for Power, Sqrt, ATan2, Pi

#include <algorithm> // for max, for_each, copy_if
#include <cmath>     // for fabs, acos
#include <cstddef>   // for size_t
#include <iostream>  // for operator<<, basic_ostream
#include <memory>    // for shared_ptr, __shared_p...

using XYZPoint = ROOT::Math::XYZPoint;

ClassImp(AtPATTERN::AtPRA);

/**
 * @brief Set initial parameters for HC.
 *
 * In track, sets GeoTheta, GeoPhi, GeoCenter, GeoRadius.
 */
void AtPATTERN::AtPRA::SetTrackInitialParameters(AtTrack &track)
{
   fTrackSeeder->SetTrackInitialParameters(track, fRadiusFitFraction, fMinHitsRadius, fMaxHitsRadius);
}

Double_t fitf(Double_t *x, Double_t *par)
{

   return par[0] + par[1] * x[0];
}

void AtPATTERN::AtPRA::OrderClustersAlongTrack(AtTrack &track)
{
   fTrackRefiner->OrderClustersAlongTrack(track);
}

void AtPATTERN::AtPRA::SelectAndMergeTracks(std::vector<AtTrack> &tracks, double vertexRadiusXY, double mergeDist,
                                            double minLabTheta)
{
   fTrackRefiner->SelectAndMergeTracks(tracks, *fTrackTransformer, fClusterRadius, fClusterDistance, vertexRadiusXY,
                                       mergeDist, minLabTheta);
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
