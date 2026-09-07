#include "AtSampleConsensus.h"
#include "AtPatternLine.h"
#include "AtPatternCircle2D.h"
#include "AtPatternTypes.h"
#include "AtPatternEvent.h"
#include "AtTrack.h"
#include <algorithm>
#include <cmath>
#include <vector>

/// ---------------------------------------------------------------------------------------------
/// AtPRA's theta determination, reproduced here so it can be SEEN.
///
/// This is AtPRA::SetTrackInitialParameters (AtPRA.cxx:36-320) replayed on the track's own hit
/// array, with every constant matching the class defaults: fRadiusFitFraction 1.0,
/// fMinHitsRadius 3, fMaxHitsRadius 1000, precluster OFF; circle RANSAC threshold 6.0 / 1000
/// iterations / minHits 0.1n; line RANSAC threshold 6.0 / minHits 0.1n / FitPattern on.
/// Nothing in the framework is modified -- the GUI just recomputes what AtPRA computed.
///
/// The point of drawing it: theta is ONE SIGN of ONE RANSAC line in this plane
/// (dir2D.X()*dir2D.Y() < 0), and the arclength is built by unwrapping the azimuth IN HIT-ARRAY
/// ORDER, so a scrambled array scrambles this plot. If the points do not lie on a line here, the
/// direction downstream is not trustworthy however good the fit looks in x-y-z.
struct PraTheta {
   std::vector<double> arc, z;   // every thetaHit, in hit-array order
   std::vector<char> inlier;     // kept by the line RANSAC?
   double dirX{0}, dirY{0};      // the fitted line direction (unit)
   double px{0}, py{0};          // a point on it
   double angleDeg{-1};          // what AtPRA would call GeoTheta
   int sign{0};
   double cx{0}, cy{0}, radius{0};
   int nUsedCircle{0}, nTotalHits{0};
   bool ok{false};
};

static PraTheta praThetaHits(const AtTrack &trk)
{
   PraTheta o;
   auto &allHits = trk.GetHitArray();
   const int nTotal = allHits.size();
   o.nTotalHits = nTotal;
   if (nTotal < 5) return o;

   // --- the hit subset AtPRA fits the circle to: the LAST nHitsForFit of the array ------------
   const double fRadiusFitFraction = 1.0; const int fMinHitsRadius = 3, fMaxHitsRadius = 1000;
   int nFromFraction = std::max(3, (int)(nTotal * fRadiusFitFraction));
   int nHitsForFit = std::max(std::min(fMinHitsRadius, nTotal),
                              std::min(nFromFraction, std::min(fMaxHitsRadius, nTotal)));
   std::vector<const AtHit *> hitsForFit;
   for (int i = std::max(0, nTotal - nHitsForFit); i < nTotal; ++i) hitsForFit.push_back(allHits[i].get());

   SampleConsensus::AtSampleConsensus rc;
   rc.SetPatternType(AtPatterns::PatternType::kCircle2D);
   rc.SetMinHitsPattern(0.1 * hitsForFit.size());
   rc.SetDistanceThreshold(6.0);
   rc.SetNumIterations(1000);
   auto circularTracks = rc.Solve(hitsForFit).GetTrackCand();
   if (circularTracks.empty()) return o;

   auto circle = std::make_unique<AtPatterns::AtPatternCircle2D>();
   circle->AtPattern::FitPattern(hitsForFit);
   const auto center = circle->GetCenter();
   o.radius = circle->GetRadius(); o.cx = center.X(); o.cy = center.Y();

   auto &hits = circularTracks.at(0).GetHitArray();     // the circle's consensus subset
   o.nUsedCircle = hits.size();
   if (hits.size() < 3) return o;

   // --- unwrap the azimuth IN HIT-ARRAY ORDER and build (arc, z) ------------------------------
   auto posPCA = hits.at(0)->GetPosition();
   auto refPosOnCircle = posPCA - center;
   const double refAng = refPosOnCircle.Phi();
   int numYCross = 0;
   auto sgn = [](double v) { return v > 0 ? 1 : (v < 0 ? -1 : 0); };
   int lastYSign = sgn(refPosOnCircle.Y());
   std::vector<AtHit> thetaHits;
   for (size_t i = 0; i < hits.size(); ++i) {
      auto pos = hits.at(i)->GetPosition();
      auto posOnCircle = pos - center;
      double angleHit = posOnCircle.Phi();
      int currYSign = sgn(posOnCircle.Y());
      if (posOnCircle.X() < 0 && lastYSign != currYSign) numYCross -= currYSign;
      lastYSign = currYSign;
      angleHit += 2 * M_PI * numYCross;
      const double a = o.radius * (refAng - angleHit);
      o.arc.push_back(a); o.z.push_back(pos.Z());
      thetaHits.emplace_back((int)i, hits.at(i)->GetPadNum(),
                             ROOT::Math::XYZPoint(a, pos.Z(), i * 1E-19), hits.at(i)->GetCharge());
   }
   o.inlier.assign(o.arc.size(), 0);

   // --- the line RANSAC whose SIGN is the whole direction decision ----------------------------
   SampleConsensus::AtSampleConsensus rt;
   rt.SetPatternType(AtPatterns::PatternType::kLine);
   rt.SetMinHitsPattern(0.1 * thetaHits.size());
   rt.SetDistanceThreshold(6.0);
   rt.SetFitPattern(true);
   auto thetaTracks = rt.Solve(thetaHits).GetTrackCand();
   if (thetaTracks.empty()) return o;
   auto *line = dynamic_cast<const AtPatterns::AtPatternLine *>(thetaTracks.at(0).GetPattern());
   if (!line) return o;
   for (auto &h : thetaTracks.at(0).GetHitArray()) {
      int id = h->GetHitID();
      if (id >= 0 && id < (int)o.inlier.size()) o.inlier[id] = 1;
   }
   auto d = line->GetDirection();
   double n = std::hypot(d.X(), d.Y());
   if (n < 1e-12) return o;
   o.dirX = d.X() / n; o.dirY = d.Y() / n;
   o.px = line->GetPoint().X(); o.py = line->GetPoint().Y();
   o.sign = (o.dirX * o.dirY < 0) ? -1 : 1;
   if (o.dirX != 0) o.angleDeg = std::acos(o.sign * std::fabs(o.dirY)) * TMath::RadToDeg();
   o.ok = true;
   return o;
}

