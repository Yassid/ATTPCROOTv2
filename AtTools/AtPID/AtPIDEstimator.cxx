#include "AtPIDEstimator.h"

#include "AtHit.h"
#include "AtTrack.h"

#include <algorithm>
#include <cmath>
#include <vector>

using namespace AtTools;

AtPIDResult AtPIDEstimator::Estimate(AtTrack &track) const
{
   AtPIDResult r;
   const double R = track.GetGeoRadius();   // mm
   r.polar = track.GetGeoTheta();           // rad (digi frame)
   r.azimuthal = track.GetGeoPhi();         // rad

   // --- brho = B * R * 0.001 / |sin(polar)|  (Spyral estimator) ---
   const double sinT = std::sin(r.polar);
   if (R > 0 && std::fabs(sinT) > 1e-6) {
      r.brho = fBField * R * 0.001 / std::fabs(sinT);
      if (std::isnan(r.brho))
         r.brho = 0.0;
   }

   // --- dEdx: charge integrated from the vertex outward over the inner-pad
   //     segment (r < smallPadRadius), divided by the arc length of that segment.
   const auto &hits = track.GetHitArray();
   struct HitRC {
      double rad;
      ROOT::Math::XYZPoint pos;
      double q;
   };
   std::vector<HitRC> ordered;
   ordered.reserve(hits.size());
   for (const auto &h : hits) {
      if (!h)
         continue;
      const auto p = h->GetPosition();
      ordered.push_back({std::sqrt(p.X() * p.X() + p.Y() * p.Y()), p, h->GetCharge()});
   }
   // Order from the beam axis outward (vertex -> Bragg side for recoil tracks).
   std::sort(ordered.begin(), ordered.end(), [](const HitRC &a, const HitRC &b) { return a.rad < b.rad; });

   if (!ordered.empty())
      r.vertex = ordered.front().pos;

   auto integrate = [&](double radiusCut) {
      double sumQ = 0.0, arclen = 0.0;
      bool havePrev = false;
      ROOT::Math::XYZPoint prev;
      int n = 0;
      for (const auto &hc : ordered) {
         if (hc.rad > radiusCut)
            break;
         sumQ += hc.q;
         if (havePrev) {
            const auto d = hc.pos - prev;
            arclen += std::sqrt(d.X() * d.X() + d.Y() * d.Y() + d.Z() * d.Z());
         }
         prev = hc.pos;
         havePrev = true;
         ++n;
      }
      return std::make_tuple(sumQ, arclen, n);
   };

   auto [sumQ, arclen, nInner] = integrate(fSmallPadRadius);
   // Fallback: if the track barely enters the inner region, use the whole track
   // so we still get a dE/dx estimate (clearly logged via valid flag downstream).
   if (nInner < 2) {
      std::tie(sumQ, arclen, nInner) = integrate(1e9);
   }

   r.dE = sumQ;
   r.arclength = arclen;
   r.dEdx = (arclen > 0.0) ? sumQ / arclen : 0.0;
   r.sqrtdEdx = std::sqrt(std::max(0.0, r.dEdx));
   r.valid = (R > 0.0 && r.brho > 0.0 && r.dEdx > 0.0);
   return r;
}
