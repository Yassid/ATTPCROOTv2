#include "AtHelix3D.h"

#include "AtHelix3DCore.h"
#include "AtHit.h"
#include "AtTrack.h"

#include <cmath>

using namespace AtTools;

AtHelix3DResult AtHelix3D::Estimate(const AtTrack &track, int spyralDirection) const
{
   AtHelix3DResult res;
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
   double brho, sq, th, path;
   if (!Helix3D::firstArcPID(r, V.end, fSmallPadRadius, brho, sq, th, path)) {
      res.failCode = 3;
      res.theta = E.theta;
      res.brho = E.brho;
      return res;
   }
   res.theta = th;
   res.brho = brho;
   res.sqrtdEdx = sq;
   res.dEdx = sq * sq;
   res.pidPath = path;
   res.valid = true;
   return res;
}
