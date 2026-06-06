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
   //     Hits are ordered by cylindrical radius (vertex -> outward for recoil
   //     tracks) and the arc length is the sum of consecutive 3D hit distances.
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
   if (nInner < 2)
      std::tie(sumQ, arclen, nInner) = integrate(1e9); // fallback: whole track

   r.dE = sumQ;
   r.arclength = arclen;

   // dEdx: truncated mean of per-step dq/ds over the inner segment. The standard
   // TPC ionization-PID estimator — dropping the high Landau tail (top fraction)
   // sharpens species bands far better than total-charge/length. Falls back to
   // total/length if too few steps.
   double radiusCut = (nInner >= 2 ? fSmallPadRadius : 1e9);
   std::vector<double> dqds;
   for (std::size_t i = 1; i < ordered.size(); ++i) {
      if (ordered[i].rad > radiusCut)
         break;
      const auto d = ordered[i].pos - ordered[i - 1].pos;
      const double ds = std::sqrt(d.X() * d.X() + d.Y() * d.Y() + d.Z() * d.Z());
      if (ds > 0.5) // ignore sub-pad steps (huge spurious dq/ds)
         dqds.push_back(ordered[i].q / ds);
   }
   if (dqds.size() >= 3) {
      std::sort(dqds.begin(), dqds.end());
      const std::size_t keep = std::max<std::size_t>(1, (std::size_t)(0.70 * dqds.size())); // lowest 70%
      double s = 0.0;
      for (std::size_t i = 0; i < keep; ++i)
         s += dqds[i];
      r.dEdx = s / keep;
   } else {
      r.dEdx = (arclen > 0.0) ? sumQ / arclen : 0.0;
   }
   r.sqrtdEdx = std::sqrt(std::max(0.0, r.dEdx));
   r.valid = (R > 0.0 && r.brho > 0.0 && r.dEdx > 0.0);
   return r;
}
