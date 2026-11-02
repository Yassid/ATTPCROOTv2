#include "AtHelix3D.h"

#include "AtHelix3DCore.h"
#include "AtHit.h"
#include "AtTrack.h"

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

using namespace AtTools;

std::vector<int> AtHelix3D::FindStrayHits(const AtTrack &track, double guardRadius) const
{
   std::vector<int> stray;
   const auto &H = track.GetHitArray();
   if ((int)H.size() < fMinHits)
      return stray;
   std::vector<Helix3D::P> h;
   h.reserve(H.size());
   for (const auto &x : H) {
      auto p = x->GetPosition();
      h.push_back({p.X(), p.Y(), p.Z(), x->GetCharge()});
   }
   Helix3D::EndOpt eo;
   eo.blobLinMin = fBlobLinMin;
   eo.blobOffCircle = fBlobOffCircle;
   auto r = Helix3D::estimate(h, fBField, fZSign, fZOffset, 5, fBinW, fMaxEdge, fMinBlob, fBlobFrac, fMaxDropFrac, eo);
   if (!r.ok)
      return stray;
   std::vector<char> kept(h.size(), 0);
   for (int i : r.order)
      kept[i] = 1;
   for (size_t i = 0; i < h.size(); ++i)
      if (!kept[i] && std::hypot(h[i].x, h[i].y) <= guardRadius)
         stray.push_back(i);
   return stray;
}

AtHelix3DResult AtHelix3D::Estimate(const AtTrack &track, int spyralDirection) const
{
   AtHelix3DResult res;
   res.trackID = track.GetTrackID();
   // AtTrack has no const accessor for the clusters; reading their count does not modify the track
   if (auto *hc = const_cast<AtTrack &>(track).GetHitClusterArray())
      res.nClusters = hc->size();
   const auto &H = track.GetHitArray();
   res.nHits = H.size();
   if (res.nHits < fMinHits) {
      res.failCode = 1;
      return res;
   }
   std::vector<Helix3D::P> h;
   h.reserve(H.size());
   for (const auto &x : H) {
      auto p = x->GetPosition();
      h.push_back({p.X(), p.Y(), p.Z(), x->GetCharge()});
   }

   Helix3D::EndOpt eo;
   eo.guardFrac = fGuardFrac;
   eo.dcaMax = fDcaMax;
   eo.rStartAxis = fRStartAxis;
   eo.blobLinMin = fBlobLinMin;
   eo.blobOffCircle = fBlobOffCircle;
   eo.winMaxFrac = fWinMaxFrac;
   eo.winMaxLen = fWinMaxLen;
   eo.hitRefit = fHitRefit;
   eo.hitRefitZ = fHitRefitZ;
   auto r = Helix3D::estimate(h, fBField, fZSign, fZOffset, 5, fBinW, fMaxEdge, fMinBlob, fBlobFrac, fMaxDropFrac, eo);
   res.nDropped = r.nDropped;
   res.length = r.length;
   res.spacing = r.nBackbone > 0 ? r.length / r.nBackbone : 0;
   if (!r.ok) {
      res.failCode = 2;
      return res;
   }

   // r.order indexes the CALLER's hits (digi frame); h is still digi frame here.
   const double zDigiA = h[r.order.front()].z, zDigiB = h[r.order.back()].z;
   int spyEnd = -1; // Spyral forward = vertex at the LOWEST z_digi
   if (spyralDirection == 0 || spyralDirection == 1)
      spyEnd = (spyralDirection == 0) == (zDigiA < zDigiB) ? 0 : 1;
   auto V = Helix3D::majorityEnd(r, spyEnd, fExpSense, fBraggThr, fSmallPadRadius, fRContained);
   res.vertexEnd = V.end;
   res.ambiguous = V.ambiguous;
   res.nFor = V.nFor;
   res.nVotes = V.nVotes;
   res.voteSense = V.sense;
   res.voteSpyral = V.spyral;
   res.voteBragg = V.bragg;
   res.braggRatio = V.braggRatio;
   const double zLabA = fZSign * zDigiA + fZOffset, zLabB = fZSign * zDigiB + fZOffset;
   res.zLabVertexEnd = V.end == 0 ? zLabA : zLabB;
   res.zLabOtherEnd = V.end == 0 ? zLabB : zLabA;
   res.backward = res.zLabVertexEnd > res.zLabOtherEnd;

   const auto &E = r.e[V.end];
   res.radius = E.R;
   res.dca = E.dca;
   res.vertex = ROOT::Math::XYZPoint(E.vx, E.vy, E.vz);
   res.guarded = E.guarded;
   res.axisConstrained = E.axisConstrained;

   const int mC = r.cen.size();
   // backbone coordinate of every ordered hit, measured from the VERTEX end
   auto sFromVertex = [&](size_t i) { return V.end == 0 ? r.s[i] : r.length - r.s[i]; };

   // --- stopping end and containment (centroids are in the internal = lab frame) ---
   if (mC >= 2) {
      const auto &cf = r.cen[V.end == 0 ? mC - 1 : 0];
      res.rStop = std::hypot(cf.x, cf.y);
      res.zLabStop = cf.z;
      res.contained = res.rStop < fContRMax && cf.z > fContZMin && cf.z < fContZMax;
   }

   // --- charge profile along the track: where does the charge peak? (macro range_dp.C, 2026-11-01) ---
   {
      const int nb = std::max(1, (int)std::ceil(r.length / 10.0));
      std::vector<double> prof(nb, 0), sm(nb, 0);
      for (size_t i = 0; i < r.order.size(); ++i)
         prof[std::min(nb - 1, std::max(0, (int)(sFromVertex(i) / 10.0)))] += h[r.order[i]].q;
      for (int b = 0; b < nb; ++b) {
         double sum = 0;
         int c = 0;
         for (int d = -1; d <= 1; ++d)
            if (b + d >= 0 && b + d < nb) {
               sum += prof[b + d];
               ++c;
            }
         sm[b] = sum / c;
      }
      const int bp = std::max_element(sm.begin(), sm.end()) - sm.begin();
      const double sPk = (bp + 0.5) * 10;
      res.peakFrac = sPk / std::max(1.0, r.length);
      double ta = 0;
      int tc = 0;
      for (int b = bp + 2; b < nb; ++b) {
         ta += sm[b];
         ++tc;
      }
      res.tailRatio = tc ? ta / tc / std::max(1e-9, sm[bp]) : 1;
      res.midPeak = res.peakFrac < fMpFpk && r.length - sPk > fMpAfter && res.tailRatio < fMpTail;
   }

   // --- arc-length polar angle over the vertex-end window (macro pid_sbpid_dp.C, 2026-11-01) ---
   if (fArcAngle && E.ok && E.nUsed >= 1 && E.R > 0) {
      const double lim = (V.end == 0 ? r.sc[E.nUsed - 1] : r.length - r.sc[mC - E.nUsed]) + 2.0;
      std::vector<std::pair<double, int>> win;
      for (size_t i = 0; i < r.order.size(); ++i)
         if (sFromVertex(i) <= lim)
            win.push_back({sFromVertex(i), r.order[i]});
      std::sort(win.begin(), win.end());
      const int en = win.size();
      if (en >= 12) {
         const int n = std::max(6, (int)(fArcFrac * en));
         double prev = std::atan2(h[win[0].second].y - E.yc, h[win[0].second].x - E.xc), ph = 0;
         double sx = 0, sy = 0, sxx = 0, sxy = 0;
         for (int i = 0; i < n && i < en; ++i) {
            const auto &hp = h[win[i].second];
            double a = std::atan2(hp.y - E.yc, hp.x - E.xc), d = a - prev;
            while (d > M_PI)
               d -= 2 * M_PI;
            while (d < -M_PI)
               d += 2 * M_PI;
            if (i)
               ph += d; // NET phase: summing |dphi| would inflate the arc with hit-order jitter
            prev = a;
            const double sl = E.R * std::abs(ph), zl = fZSign * hp.z + fZOffset;
            sx += sl;
            sy += zl;
            sxx += sl * sl;
            sxy += sl * zl;
         }
         const int m = std::min(n, en);
         const double den = m * sxx - sx * sx;
         if (std::abs(den) > 1e-9) {
            res.thetaArc = std::atan2(1.0, (m * sxy - sx * sy) / den);
            res.brhoArc = fBField * E.R / 1000.0 / std::sin(res.thetaArc);
         }
      }
   }
   double brho, sq, th, path;
   if (!Helix3D::firstArcPID(r, V.end, fSmallPadRadius, brho, sq, th, path)) {
      res.failCode = 3;
      res.theta = E.theta;
      res.brho = E.brho;
      return res;
   }
   res.theta = th;
   res.brho = brho;
   res.sqrtdEdx = sq * fDedxGain;
   res.dEdx = res.sqrtdEdx * res.sqrtdEdx;
   res.pidPath = path;
   res.valid = true;
   return res;
}
