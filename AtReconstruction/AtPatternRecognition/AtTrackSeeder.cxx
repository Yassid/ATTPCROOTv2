#include "AtTrackSeeder.h"

#include "AtHit.h"
#include "AtPatternCircle2D.h"
#include "AtPatternEvent.h"
#include "AtPatternLine.h"
#include "AtPatternTypes.h"
#include "AtSampleConsensus.h"
#include "AtTrack.h"

#include <FairLogger.h>

#include <Math/Point3D.h>
#include <Math/Vector2D.h>
#include <TGraph.h>
#include <TMath.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <memory>
#include <vector>

void AtPATTERN::AtTrackSeeder::SetTrackInitialParameters(AtTrack &track, double radiusFitFraction, int minHitsRadius,
                                                         int maxHitsRadius)
{
   auto &allHits = track.GetHitArray();
   int nTotal = allHits.size();
   int nFromFraction = std::max(3, static_cast<int>(nTotal * radiusFitFraction));
   int nHitsForFit =
      std::max(std::min(minHitsRadius, nTotal), std::min(nFromFraction, std::min(maxHitsRadius, nTotal)));

   std::vector<const AtHit *> hitsForFit;
   int startIdx = std::max(0, nTotal - nHitsForFit);
   for (int i = startIdx; i < nTotal; i++)
      hitsForFit.push_back(allHits[i].get());

   LOG(debug) << "AtTrackSeeder: circle fit using " << hitsForFit.size() << "/" << nTotal << " hits";

   SampleConsensus::AtSampleConsensus ransacSmoothRadius;
   ransacSmoothRadius.SetPatternType(AtPatterns::PatternType::kCircle2D);
   ransacSmoothRadius.SetMinHitsPattern(0.1 * hitsForFit.size());
   ransacSmoothRadius.SetDistanceThreshold(6.0);
   ransacSmoothRadius.SetNumIterations(1000);
   auto circularTracks = ransacSmoothRadius.Solve(hitsForFit).GetTrackCand();

   if (circularTracks.empty())
      return;

   auto &hits = circularTracks.at(0).GetHitArray();

   auto circle = std::make_unique<AtPatterns::AtPatternCircle2D>();
   circle->AtPattern::FitPattern(hitsForFit);

   auto center = circle->GetCenter();
   auto radius = circle->GetRadius();

   track.SetGeoCenter({center.X(), center.Y()});
   track.SetGeoRadius(radius);
   track.SetPattern(circle->Clone());

   circularTracks.at(0).SetPattern(std::move(circle));

   std::vector<double> whit;
   std::vector<double> arclength;

   auto arclengthGraph = std::make_unique<TGraph>();

   auto posPCA = hits.at(0)->GetPosition();
   auto refPosOnCircle = posPCA - center;
   auto refAng = refPosOnCircle.Phi();

   std::vector<AtHit> thetaHits;

   int numYCross = 0;
   int lastYSign = (0 < refPosOnCircle.Y()) - (refPosOnCircle.Y() < 0);

   for (size_t i = 0; i < hits.size(); ++i) {
      auto pos = hits.at(i)->GetPosition();
      auto posOnCircle = pos - center;
      auto angleHit = posOnCircle.Phi();

      int currYSign = (0 < posOnCircle.Y()) - (posOnCircle.Y() < 0);
      if (posOnCircle.X() < 0 && lastYSign != currYSign)
         numYCross -= currYSign;
      lastYSign = currYSign;
      angleHit += 2 * M_PI * numYCross;

      whit.push_back(angleHit);
      arclength.push_back((radius * (refAng - whit.at(i))));

      arclengthGraph->SetPoint(arclengthGraph->GetN(), arclength.at(i), pos.Z());

      Double_t xPos = arclength.at(i);
      Double_t yPos = pos.Z();
      Double_t zPos = i * 1E-19;
      thetaHits.emplace_back(i, hits.at(i)->GetPadNum(), ROOT::Math::XYZPoint(xPos, yPos, zPos), hits.at(i)->GetCharge());
   }

   Double_t angle = 0.0;
   Double_t phi0 = 0.0;

   try {
      if (!thetaHits.empty()) {
         std::vector<AtTrack> thetaTracks;

         SampleConsensus::AtSampleConsensus ransacTheta;
         ransacTheta.SetPatternType(AtPatterns::PatternType::kLine);
         ransacTheta.SetMinHitsPattern(0.1 * thetaHits.size());
         ransacTheta.SetDistanceThreshold(6.0);
         ransacTheta.SetFitPattern(true);
         thetaTracks = ransacTheta.Solve(thetaHits).GetTrackCand();

         if (!thetaTracks.empty()) {
            auto line = dynamic_cast<const AtPatterns::AtPatternLine *>(thetaTracks.at(0).GetPattern());
            LOG(info) << "Track ID: " << track.GetTrackID() << " with " << track.GetHitArray().size() << " hits. Fit "
                      << thetaTracks[0].GetHitArray().size() << "/" << thetaHits.size() << " hits in first of "
                      << thetaTracks.size() << " tracks";

            auto dirTheta = line->GetDirection();
            auto dir2D = ROOT::Math::XYVector(dirTheta.X(), dirTheta.Y()).Unit();

            int sign = dir2D.X() * dir2D.Y() < 0 ? -1 : 1;

            if (dir2D.X() != 0)
               angle = acos(sign * fabs(dir2D.Y())) * TMath::RadToDeg();

            LOG(info) << "Setting theta geo to: " << angle << " with ransac direction: " << dir2D;
            for (int i = 1; i < thetaTracks.size(); ++i) {
               auto ranLine = dynamic_cast<const AtPatterns::AtPatternLine *>(thetaTracks.at(i).GetPattern());
               auto ranDir = ROOT::Math::XYVector(ranLine->GetDirection().X(), ranLine->GetDirection().Y()).Unit();
               LOG(info) << "RANSAC direction " << i << " with " << thetaTracks[i].GetHitArray().size() << "/"
                         << thetaHits.size() << ": " << ranDir;
            }

            auto temp2 = posPCA - center;
            phi0 = TMath::ATan2(temp2.Y(), temp2.X());
         }
         track.SetGeoTheta(angle * TMath::Pi() / 180.0);
         track.SetGeoPhi(phi0);
      } else {
         LOG(info) << "Track ID: " << track.GetTrackID() << " with " << track.GetHitArray().size()
                   << " hits does not have enough theta hits to get angle";
      }
   } catch (std::exception &e) {
      std::cout << " AtTrackSeeder::SetTrackInitialParameters - Exception caught : " << e.what() << "\n";
   }
}
