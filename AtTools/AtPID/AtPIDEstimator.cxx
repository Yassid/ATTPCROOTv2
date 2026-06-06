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

   // --- dEdx: Spyral recipe = total charge integrated over the inner-pad
   //     segment (r < smallPadRadius) / arc length of that segment. Computed on
   //     the PRA CLUSTERS (charge-weighted, denoised) rather than raw hits — this
   //     is what makes the bands sharp. Clusters carry the summed charge in
   //     GetCharge(); we order them by cylindrical radius (vertex -> outward) and
   //     sum consecutive 3D distances for the arc length. Falls back to raw hits
   //     if the track has no cluster array.
   struct PtRC {
      double rad;
      ROOT::Math::XYZPoint pos;
      double q;
   };
   std::vector<PtRC> pts;
   const auto *clusters = track.GetHitClusterArray();
   if (clusters && clusters->size() >= 3) {
      pts.reserve(clusters->size());
      for (const auto &c : *clusters) {
         auto p = c.GetPositionCharge(); // charge-weighted centroid
         if (p.X() == 0 && p.Y() == 0 && p.Z() == 0)
            p = c.GetPosition();
         pts.push_back({std::sqrt(p.X() * p.X() + p.Y() * p.Y()), ROOT::Math::XYZPoint(p.X(), p.Y(), p.Z()),
                        c.GetCharge()});
      }
   } else {
      const auto &hits = track.GetHitArray();
      pts.reserve(hits.size());
      for (const auto &h : hits) {
         if (!h)
            continue;
         const auto p = h->GetPosition();
         pts.push_back({std::sqrt(p.X() * p.X() + p.Y() * p.Y()), p, h->GetCharge()});
      }
   }
   std::sort(pts.begin(), pts.end(), [](const PtRC &a, const PtRC &b) { return a.rad < b.rad; });
   if (!pts.empty())
      r.vertex = pts.front().pos;
   r.nClusters = static_cast<int>(pts.size());

   auto integrate = [&](double radiusCut) {
      double sumQ = 0.0, arclen = 0.0;
      bool havePrev = false;
      ROOT::Math::XYZPoint prev;
      int n = 0;
      for (const auto &pt : pts) {
         if (pt.rad > radiusCut)
            break;
         sumQ += pt.q;
         if (havePrev) {
            const auto d = pt.pos - prev;
            arclen += std::sqrt(d.X() * d.X() + d.Y() * d.Y() + d.Z() * d.Z());
         }
         prev = pt.pos;
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
   r.dEdx = (arclen > 0.0) ? sumQ / arclen : 0.0;
   r.sqrtdEdx = std::sqrt(std::max(0.0, r.dEdx));

   // --- Native pid_dev.C variables: mean cluster charge over the first 80% of
   //     clusters from the Bragg end (forward tracks: high-index end; backward:
   //     low-index end), in the cluster array's intrinsic (Z-sorted) order. This
   //     is length-division-free (no noisy arclength) and emphasizes the
   //     Bragg-rise region where species separate best. elossTrunc drops the
   //     high Landau tail.
   if (clusters && clusters->size() >= 3) {
      const double thetaDeg = r.polar * 180.0 / M_PI;
      const int nc = static_cast<int>(clusters->size());
      const int maxN = std::max(1, static_cast<int>(0.8 * nc));
      std::vector<double> q;
      q.reserve(maxN);
      if (thetaDeg < 90.0) {
         for (int k = nc - 1; k >= 0 && static_cast<int>(q.size()) < maxN; --k)
            q.push_back((*clusters)[k].GetCharge());
      } else {
         for (int k = 0; k < nc && static_cast<int>(q.size()) < maxN; ++k)
            q.push_back((*clusters)[k].GetCharge());
      }
      double s = 0.0;
      for (double v : q)
         s += v;
      r.elossMean = q.empty() ? 0.0 : s / q.size();
      std::sort(q.begin(), q.end());
      const int keep = std::max(1, static_cast<int>(0.7 * q.size()));
      double st = 0.0;
      for (int k = 0; k < keep; ++k)
         st += q[k];
      r.elossTrunc = st / keep;
   }

   r.valid = (R > 0.0 && r.brho > 0.0 && r.dEdx > 0.0);
   return r;
}
