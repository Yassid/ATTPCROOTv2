#include "AtGenfitter.h"

#include "AtFitTrackMetadata.h"
#include "AtFittedTrack.h"
#include "AtHitCluster.h"
#include "AtParticleID.h"
#include "AtPatternEvent.h"
#include "AtSpacePointMeasurement.h"
#include "AtSpyralPID.h"
#include "AtTrack.h"
#include "AtTrackingEvent.h"

#include <FairLogger.h>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <TClonesArray.h>
#include <TDatabasePDG.h>
#include <TGeoManager.h>
#include <TMatrixDSym.h>
#include <TVector3.h>

#include "ConstField.h"
#include "Exception.h"
#include "FieldManager.h"
#include "FitStatus.h"
#include "KalmanFitterRefTrack.h"
#include "MaterialEffects.h"
#include "MeasuredStateOnPlane.h"
#include "MeasurementFactory.h"
#include "MeasurementProducer.h"
#include "RKTrackRep.h"
#include "TGeoMaterialInterface.h"
#include "Track.h"
#include "TrackCand.h"

#include <algorithm>
#include <cmath>

ClassImp(EventFit::AtGenfitter);

using XYZPoint = ROOT::Math::XYZPoint;
using XYZVector = ROOT::Math::XYZVector;

namespace EventFit {

AtGenfitter::AtGenfitter(Double_t bFieldTesla, Int_t pdg, Double_t massAmu, Int_t Z, std::string eLossFile,
                         Bool_t noMatEffects, Int_t minIter, Int_t maxIter)
   : fBField(bFieldTesla), fPDG(pdg), fMassAmu(massAmu), fZ(Z), fELossFile(std::move(eLossFile)),
     fNoMatEffects(noMatEffects), fMinIter(minIter), fMaxIter(maxIter)
{
   fHitClusterArray = new TClonesArray("AtHitCluster");
   fGenfitTrackArray = new TClonesArray("genfit::Track");
   fMeasurementProducer = new genfit::MeasurementProducer<AtHitCluster, genfit::AtSpacepointMeasurement>(fHitClusterArray);
   fMeasurementFactory = new genfit::MeasurementFactory<genfit::AbsMeasurement>();
}

AtGenfitter::~AtGenfitter()
{
   delete fGenfitTrackArray;
   delete fHitClusterArray;
   // NOTE: do NOT delete fMeasurementProducer here. Init() registers it with the
   // MeasurementFactory via addProducer(), and MeasurementFactory::clear() (called by
   // ~MeasurementFactory) deletes every registered producer. Deleting it here too was a
   // double-free that crashed when more than one AtGenfitter was constructed/destroyed.
   delete fMeasurementFactory;
}

void AtGenfitter::Init()
{
   if (fInit)
      return;
   fKalmanFitter = std::make_shared<genfit::KalmanFitterRefTrack>();
   fKalmanFitter->setMinIterations(fMinIter);
   fKalmanFitter->setMaxIterations(fMaxIter);
   fMeasurementFactory->addProducer(fTPCDetID, fMeasurementProducer);

   // signed B along z, in kGauss (genfit convention used by AtGenfit: 10 * Tesla)
   genfit::FieldManager::getInstance()->init(new genfit::ConstField(0., 0., 10.0 * fBField));
   genfit::MaterialEffects *mat = genfit::MaterialEffects::getInstance();
   mat->setEnergyLossBrems(false);
   mat->setNoiseBrems(false);
   mat->setNoEffects(fNoMatEffects);
   mat->init(new genfit::TGeoMaterialInterface());

   // genfit's default (true) makes MaterialEffects::stepper() walk forward across every
   // boundary that separates volumes of the SAME material, merging them into one step. That
   // is an optimization for sparse trackers whose material is split over many volumes. The
   // AT-TPC active volume is a SINGLE homogeneous gas tube, so the "material after" test can
   // never fire and the loop simply runs to its 100-iteration cap, calling findNextBoundary
   // (itself up to 300 RKPropagate-from-origin iterations plus TGeo navigation) ~83 times per
   // RK step for a step limit that is not even the binding one. Measured on a1954 run_0055:
   // 113.1M findNextBoundary calls -> 1.37M, 196 s -> 19.5 s CPU per 300 events (10x), with
   // bit-identical output (max |dKE| 4.1e-7 MeV, max |dTheta| 2.0e-6 deg, same track count,
   // same genfit exception counts). Turning it off also stops the OOM kills on long runs.
   mat->ignoreBoundariesBetweenEqualMaterials(false);

   // Energy-loss mode. genfit's parameterization-only mode (useEnergyLossParam) is ONLY
   // valid when a dE/dx curve has been loaded: MaterialEffects::maxKinEnergy_ starts at 0
   // and is set solely by setEnergyLossFile(), and dEdx() throws a FATAL exception for
   // E_kin > maxKinEnergy_. Calling it without a table therefore threw on every single
   // step, which silently truncated genfit's reference track after ~2 measurements
   // (ndf <= 0) for ~87 % of tracks -- the "material effects are unusable" instability.
   // Without a table we leave genfit in Param+Bethe-Bloch mode, which is valid down to
   // beta*gamma = 0.05 (~1.2 MeV protons) and falls back to a 0 parameterization term
   // below that.
   if (!fELossFile.empty()) {
      mat->setEnergyLossFile(fELossFile, fPDG);
      mat->useEnergyLossParam();
   }

   // PDG defs (genfit needs ion entries)
   TDatabasePDG *db = TDatabasePDG::Instance();
   if (!db->GetParticle(1000010020))
      db->AddParticle("Deuteron", "Deuteron", 1.875612928, kTRUE, 0, 3, "Ion", 1000010020);
   if (!db->GetParticle(1000010030))
      db->AddParticle("Triton", "Triton", 2.80943211, kFALSE, 0, 3, "Ion", 1000010030);

   // optional upstream PID gate: compute Spyral PID per track, fit only gated species
   if (fUsePIDGate) {
      fSpyralPID = std::make_unique<AtTools::AtSpyralPID>();
      fSpyralPID->SetBField(std::abs(fBField)); // brho uses |B|, matches AtPIDTask
      auto loaded = AtTools::AtParticleID::LoadJSON(fPidGateFile);
      if (loaded.GetCut().IsValid()) {
         fPidGate = std::make_unique<AtTools::AtParticleID>(std::move(loaded));
         LOG(info) << "AtGenfitter PID gate: '" << fPidGate->GetName() << "' from " << fPidGateFile;
      } else {
         LOG(error) << "AtGenfitter: could not load PID gate '" << fPidGateFile << "' -> gating DISABLED";
         fUsePIDGate = kFALSE;
      }
   }

   fInit = kTRUE;
   LOG(info) << "AtGenfitter init: B=" << fBField << " T, pdg=" << fPDG << ", iter " << fMinIter << "-" << fMaxIter
             << ", PIDgate=" << (fUsePIDGate ? "on" : "off");
}

AtFittedTrack *AtGenfitter::GetFittedTrack(AtTrack *track, AtFitMetadata * /*fitMetadata*/, AtRawEvent * /*rawEvent*/,
                                           AtEvent * /*event*/)
{
   if (!fInit) // AtFitterTask does not call Init(); self-init on first use (geometry is loaded by now)
      Init();
   auto *hc = track->GetHitClusterArray();
   if (hc == nullptr || (int)hc->size() < fMinClusters)
      return nullptr;
   const int n = hc->size();

   // upstream PID gate: skip non-gated species before doing the fit
   if (fUsePIDGate && fPidGate) {
      auto pidRes = fSpyralPID->Estimate(*track);
      if (!pidRes.valid || !fPidGate->IsInside(pidRes.sqrtdEdx, pidRes.brho))
         return nullptr;
   }

   // --- lab-frame positions: z_lab = ZPadPlane - z_digi (uniform), x,y as-is ---
   std::vector<XYZPoint> pos(n);
   for (int i = 0; i < n; ++i) {
      auto p = hc->at(i).GetPosition();
      pos[i] = XYZPoint(p.X(), p.Y(), fZPadPlane - p.Z());
   }

   // --- order clusters by the drift coordinate z (the natural AT-TPC ordering) ---
   // z comes from drift time and is monotonic ALONG the trajectory even for curling
   // /spiraling tracks, so sorting by z sequences multi-turn tracks correctly (where
   // the old vertex->outward index heuristic scrambled them and the fit diverged).
   std::vector<int> order(n);
   for (int i = 0; i < n; ++i) order[i] = i;
   std::sort(order.begin(), order.end(), [&](int a, int b) { return pos[a].Z() < pos[b].Z(); });

   // vertex = highest z_digi end (= lowest z_lab, i.e. order.front() after the ascending
   // z sort) -- the deterministic AT-TPC convention used by PRA/the validated UKF. Using
   // an axis-distance heuristic instead flips near-axis forward tracks (both ends hug the
   // beam), so we do NOT guess the end; z ordering already put the vertex first.
   double dVtx = 1e18;
   for (int i = 0; i < n; ++i)
      dVtx = std::min(dVtx, std::hypot(pos[i].X(), pos[i].Y()));
   if (dVtx > fVertexAxisMaxDist)
      return nullptr; // never comes near the beam axis -> not a beam-vertex track

   // Vertex end: default = lowest-z_lab end (order.front()), seed dir -> +z_lab. Correct
   // for forward tracks but REVERSES backward ones (genfit reflects them into the forward
   // hemisphere). With fBackwardSeedFix, tracks PRA-tagged backward (GeoTheta > 90 deg)
   // are seeded from the highest-z_lab end with the momentum pointing INTO the track
   // (-z_lab). Forward tracks are byte-for-byte unchanged (backwardSeed stays false).
   bool backwardSeed = fBackwardSeedFix && (track->GetGeoTheta() > M_PI / 2.0);
   int iVtx = backwardSeed ? order.back() : order.front();
   int iSec = (order.size() > 1) ? (backwardSeed ? order[order.size() - 2] : order[1]) : iVtx;

   XYZPoint vPos = pos[iVtx];
   XYZPoint sPos = pos[iSec];
   XYZVector dir = (sPos - vPos);
   if (dir.R() < 1e-6)
      return nullptr;
   double theta = dir.Theta();
   double phi = dir.Phi();

   // momentum magnitude from Brho (circle radius from PRA)
   double radius_m = track->GetGeoRadius() / 1000.0; // mm -> m
   double sinth = std::sin(theta);
   if (std::abs(sinth) < 1e-3) sinth = (sinth < 0 ? -1e-3 : 1e-3);
   double brho = std::abs(fBField) * radius_m / std::abs(sinth);          // T*m
   double p_GeV = 0.299792458 * brho;                                      // |q|=e
   if (!(p_GeV > 0) || p_GeV > 100)
      return nullptr;

   // --- build the genfit track from the (ordered, lab-frame) clusters ---
   // Per-cluster measurement covariance (mm^2; AtSpacepointMeasurement -> cm^2): a
   // diffusion model where the transverse (x,y) and drift-z variance grow with the
   // drift distance L (sigma^2 = sigma0^2 + D^2 * L_cm), optionally scaled by cluster
   // charge. With the default diffusion coefficients (0) and fZLongFactor=2 this
   // reduces EXACTLY to the previous flat (s2, s2, 2*s2) covariance.
   const double s2 = fMeasSigmaMM * fMeasSigmaMM;
   genfit::TrackCand trackCand;
   fHitClusterArray->Clear("C");
   for (int oi = 0; oi < n; ++oi) {
      const int ci = order[oi];
      double Ldrift_cm = std::max(0.0, (fZPadPlane - pos[ci].Z()) / 10.0);   // drift distance (cm)
      double varT = s2 + fDiffTransMM * fDiffTransMM * Ldrift_cm;            // transverse (mm^2)
      double varZ = fZLongFactor * s2 + fDiffLongMM * fDiffLongMM * Ldrift_cm; // drift-z (mm^2)
      if (fChargeRefForCov > 0) {
         double q = hc->at(ci).GetCharge();
         double scale = (q > 0) ? std::clamp(fChargeRefForCov / q, 0.25, 4.0) : 1.0;
         varT *= scale;
         varZ *= scale;
      }
      TMatrixDSym measCov(3);
      measCov.Zero();
      measCov(0, 0) = varT;
      measCov(1, 1) = varT;
      measCov(2, 2) = varZ;

      AtHitCluster cl = hc->at(ci);
      cl.SetPosition({pos[ci].X(), pos[ci].Y(), pos[ci].Z()}); // lab-frame
      cl.SetCovMatrix(measCov);
      int idx = fHitClusterArray->GetEntriesFast();
      new ((*fHitClusterArray)[idx]) AtHitCluster(cl);
      trackCand.addHit(fTPCDetID, idx);
   }

   TVector3 posSeed(vPos.X() / 10.0, vPos.Y() / 10.0, vPos.Z() / 10.0); // cm
   TVector3 momSeed;
   momSeed.SetMagThetaPhi(p_GeV, theta, phi);

   TMatrixDSym covSeed(6);
   covSeed.Zero();
   for (int d = 0; d < 3; ++d) covSeed(d, d) = 0.01;        // 1 mm^2 in cm^2
   for (int d = 3; d < 6; ++d) covSeed(d, d) = std::pow(0.3 * p_GeV, 2) + 1e-6;
   trackCand.setCovSeed(covSeed);
   trackCand.setPosMomSeed(posSeed, momSeed, fZ);
   trackCand.setPdgCode(fPDG);

   TVector3 posRes, momRes;
   TMatrixDSym covRes(6);
   double chi2 = -1, ndf = -1;
   bool converged = false;
   std::vector<XYZPoint> fittedPts; // genfit fitted position at each cluster (lab-frame mm)

   // build a fresh genfit track from trackCand, fit it, and extract results; returns
   // false on any genfit failure (the RK loops THROW on non-convergence rather than
   // spinning, so this can fail but never hang).
   auto doFit = [&]() -> bool {
      // Free the previous fit's genfit::Track (its track points, reps and fitted states).
      // Without this the array grows by one heavy Track per fit and is never released ->
      // multi-GB growth over a full run that OOM-crashes the VM (esp. multi-parallel).
      // Results are copied into locals below, so nothing downstream references it.
      fGenfitTrackArray->Delete();
      auto *gfTrack = new ((*fGenfitTrackArray)[fGenfitTrackArray->GetEntriesFast()])
                         genfit::Track(trackCand, *fMeasurementFactory);
      gfTrack->addTrackRep(new genfit::RKTrackRep(fPDG));
      auto *rep = dynamic_cast<genfit::RKTrackRep *>(gfTrack->getTrackRep(0));
      try {
         fKalmanFitter->processTrackWithRep(gfTrack, rep, false);
         auto *st = gfTrack->getFitStatus(rep);
         converged = st->isFitConverged();
         chi2 = st->getChi2();
         ndf = st->getNdf();
         genfit::MeasuredStateOnPlane fitState = gfTrack->getFittedState();
         fitState.getPosMomCov(posRes, momRes, covRes);
         fittedPts.clear();
         int npts = gfTrack->getNumPointsWithMeasurement();
         for (int ip = 0; ip < npts; ++ip) {
            try {
               genfit::MeasuredStateOnPlane mop = gfTrack->getFittedState(ip, rep);
               TVector3 pp = mop.getPos(); // cm, lab frame
               fittedPts.emplace_back(pp.X() * 10.0, pp.Y() * 10.0, pp.Z() * 10.0); // mm
            } catch (...) {
            }
         }
         return true;
      } catch (...) {
         return false;
      }
   };

   // Material-effects provenance of the result we end up keeping. Starts as "whatever was
   // requested"; the fallback below flips it if the kept fit came from the no-material retry.
   bool usedMatEffects = !fNoMatEffects;
   bool matFallback = false;

   bool ok = doFit();
   if (!ok && !fNoMatEffects && fMatEffectsFallback) {
      // material-effects fit failed (e.g. a stopping multi-turn spiral whose RK
      // extrapolation throws). Retry this track WITHOUT material effects so it still
      // gets a (flagged) fit instead of vanishing. Clean tracks keep the better fit.
      genfit::MaterialEffects *me = genfit::MaterialEffects::getInstance();
      me->setNoEffects(true);
      ok = doFit();
      me->setNoEffects(false); // restore for the next (clean) track
      if (ok) {
         // This track is NOT from the same model as its neighbours. It must be flagged, or
         // the two populations blend into one spectrum and inflate the width.
         usedMatEffects = false;
         matFallback = true;
      }
   }
   if (!ok) {
      LOG(debug) << "AtGenfitter: fit failed on track " << track->GetTrackID();
      return nullptr;
   }

   double p = momRes.Mag() * 1000.0;          // GeV -> MeV
   double mass = fMassAmu * 931.49410242;     // MeV
   if (!std::isfinite(p) || p <= 0 || p > 1e6)
      return nullptr;
   double KE = std::sqrt(p * p + mass * mass) - mass;

   // drop unphysical near-beam / backward tracks by fitted polar angle
   double thetaFitDeg = momRes.Theta() * 180.0 / M_PI;
   if (thetaFitDeg < fThetaMinDeg || thetaFitDeg > fThetaMaxDeg)
      return nullptr;

   auto owner = std::make_unique<AtFittedTrack>();
   AtFittedTrack *ft = owner.get();
   ft->SetTrackID(track->GetTrackID());
   ft->SetKinematics(KE, momRes.Theta(), momRes.Phi());
   ft->SetKinematicsXtr(KE, momRes.Theta(), momRes.Phi());
   ft->SetVertex(XYZVector(posRes.X() * 10.0, posRes.Y() * 10.0, posRes.Z() * 10.0)); // cm -> mm
   ft->SetParticleInfo(std::to_string(fZ), fZ, fMassAmu);
   ft->SetSmoothedPositions(std::move(fittedPts)); // fitted trajectory for display/QA

   auto meta = std::make_unique<AtFitTrackMetadata>();
   meta->SetTrackID(track->GetTrackID());
   meta->SetChi2(chi2);
   meta->SetNdf(ndf);
   double chi2ndf = (ndf > 0) ? chi2 / ndf : 1e9;
   meta->SetGoodFit(converged && ndf > 0 && chi2ndf < fChi2NdfMax && KE > fKEMin && KE < fKEMax);
   meta->SetPValue(0);
   meta->SetFitConverged(converged);
   meta->SetMatEffects(usedMatEffects);
   meta->SetMatEffectsFallback(matFallback);
   ft->SetTrackMetadata(std::move(meta));

   return owner.release();
}

// ─── Continuity merging ──────────────────────────────────────────────────────────
namespace {
// 3D distance (mm) between two cluster positions.
double clusterDist(const AtHitCluster &a, const AtHitCluster &b)
{
   auto pa = a.GetPosition();
   auto pb = b.GetPosition();
   return std::sqrt(std::pow(pa.X() - pb.X(), 2) + std::pow(pa.Y() - pb.Y(), 2) + std::pow(pa.Z() - pb.Z(), 2));
}
// A track's two ends classified by distance to the beam axis (XY radius): the
// vertex-end (closest to the axis) and the far-end (farthest). This is the
// physically meaningful split for continuity: PRA fragments one track into pieces
// joined far-end -> next piece, whereas two distinct tracks sharing the production
// vertex meet vertex-end <-> vertex-end (and diverge). Robust to unordered clusters.
struct TrackEnds {
   int vtx, far;
};
TrackEnds axisEnds(const std::vector<AtHitCluster> &cl)
{
   int iv = 0, ifar = 0;
   double rmin = 1e18, rmax = -1.0;
   for (int i = 0; i < (int)cl.size(); ++i) {
      auto p = cl[i].GetPosition();
      double r = std::hypot(p.X(), p.Y());
      if (r < rmin) { rmin = r; iv = i; }
      if (r > rmax) { rmax = r; ifar = i; }
   }
   return {iv, ifar};
}
} // namespace

std::vector<AtTrack> AtGenfitter::MergeContinuousTracks(std::vector<AtTrack> tracks) const
{
   if (!fMergeContinuity || tracks.size() < 2)
      return tracks;

   // Two fragments merge when their PRA circles are compatible (centre + radius) AND
   // their nearest endpoints are within fMergeGapMM in 3D.
   auto shouldMerge = [&](AtTrack &A, AtTrack &B) -> bool {
      auto *clA = A.GetHitClusterArray();
      auto *clB = B.GetHitClusterArray();
      if (clA->empty() || clB->empty())
         return false;
      auto cA = A.GetGeoCenter();
      auto cB = B.GetGeoCenter();
      if (std::hypot(cA.first - cB.first, cA.second - cB.second) > fMergeCenterDist)
         return false;
      double rA = A.GetGeoRadius(), rB = B.GetGeoRadius();
      double rMax = std::max(rA, rB);
      if (rMax > 0 && std::abs(rA - rB) > fMergeRadiusFrac * rMax)
         return false;

      // Continuity is a far-end -> next-piece junction. The vertex-end <-> vertex-end
      // junction is two distinct tracks sharing the production vertex (they diverge),
      // NOT a fragmented track, so it must NOT merge. Take the smallest gap over the
      // three ALLOWED junctions (farA-vtxB, farA-farB, vtxA-farB); the vtxA-vtxB pair
      // is excluded. Merge only if a far-end participates and the gap is small.
      TrackEnds eA = axisEnds(*clA);
      TrackEnds eB = axisEnds(*clB);
      double gap = std::min({clusterDist((*clA)[eA.far], (*clB)[eB.vtx]), clusterDist((*clA)[eA.far], (*clB)[eB.far]),
                             clusterDist((*clA)[eA.vtx], (*clB)[eB.far])});
      double vtxGap = clusterDist((*clA)[eA.vtx], (*clB)[eB.vtx]);
      if (vtxGap < gap) // the closest approach is vertex-to-vertex -> shared vertex, diverging tracks
         return false;
      return gap < fMergeGapMM;
   };

   std::vector<bool> used(tracks.size(), false);
   std::vector<AtTrack> out;
   for (size_t i = 0; i < tracks.size(); ++i) {
      if (used[i])
         continue;
      AtTrack merged = tracks[i];
      used[i] = true;
      bool grew = true;
      while (grew) { // grow the seed track until no more fragments attach (transitive chains)
         grew = false;
         for (size_t j = 0; j < tracks.size(); ++j) {
            if (used[j] || !shouldMerge(merged, tracks[j]))
               continue;
            auto *dst = merged.GetHitClusterArray();
            for (const auto &c : *tracks[j].GetHitClusterArray())
               dst->push_back(c); // concatenate clusters; the fit's drift-z sort reorders
            used[j] = true;
            grew = true;
         }
      }
      out.push_back(std::move(merged));
   }
   if (out.size() != tracks.size())
      LOG(debug) << "AtGenfitter continuity merge: " << tracks.size() << " -> " << out.size() << " tracks";
   return out;
}

void AtGenfitter::FitEvent(AtTrackingEvent *trackingEvent, AtPatternEvent *patternEvent, AtFitMetadata *fitMetadata,
                           AtRawEvent *rawEvent, AtEvent *event)
{
   if (!fMergeContinuity) { // default path: unchanged base behaviour
      AtFitter::FitEvent(trackingEvent, patternEvent, fitMetadata, rawEvent, event);
      return;
   }
   if (trackingEvent == nullptr || patternEvent == nullptr) {
      LOG(error) << "AtGenfitter::FitEvent: null tracking/pattern event";
      return;
   }
   // Merge a COPY of the pattern tracks (the shared AtPatternEvent stays intact for any
   // parallel fitter task on the same event), then run the normal per-track fit loop.
   std::vector<AtTrack> tracks = MergeContinuousTracks(patternEvent->GetTrackCand());
   if (tracks.empty())
      return;
   trackingEvent->SetTrackArray(&tracks);
   for (auto track : tracks) {
      std::unique_ptr<AtFittedTrack> fittedTrack(GetFittedTrack(&track, fitMetadata, rawEvent, event));
      if (fittedTrack)
         trackingEvent->AddFittedTrack(std::move(fittedTrack));
   }
}

} // namespace EventFit
