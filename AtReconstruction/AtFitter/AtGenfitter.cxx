#include "AtGenfitter.h"

#include "AtFitTrackMetadata.h"
#include "AtFittedTrack.h"
#include "AtHitCluster.h"
#include "AtParticleID.h"
#include "AtSpacePointMeasurement.h"
#include "AtSpyralPID.h"
#include "AtTrack.h"

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
   delete fMeasurementProducer;
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
   mat->useEnergyLossParam();
   mat->init(new genfit::TGeoMaterialInterface());
   if (!fELossFile.empty())
      mat->setEnergyLossFile(fELossFile, fPDG);

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

   XYZPoint vPos = pos[order.front()];
   XYZPoint sPos = pos[order.size() > 1 ? order[1] : order.front()];
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

   // measurement covariance fed to genfit (mm^2; AtSpacepointMeasurement converts to cm^2)
   TMatrixDSym measCov(3);
   measCov.Zero();
   double s2 = fMeasSigmaMM * fMeasSigmaMM;
   measCov(0, 0) = s2;
   measCov(1, 1) = s2;
   measCov(2, 2) = 2.0 * s2; // z (drift) a bit looser

   // --- build the genfit track from the (ordered, lab-frame) clusters ---
   genfit::TrackCand trackCand;
   fHitClusterArray->Clear("C");
   for (int oi = 0; oi < n; ++oi) {
      AtHitCluster cl = hc->at(order[oi]);
      cl.SetPosition({pos[order[oi]].X(), pos[order[oi]].Y(), pos[order[oi]].Z()}); // lab-frame
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

   bool ok = doFit();
   if (!ok && !fNoMatEffects) {
      // material-effects fit failed (e.g. a stopping multi-turn spiral whose RK
      // extrapolation throws). Retry this track WITHOUT material effects so it still
      // gets a (flagged) fit instead of vanishing. Clean tracks keep the better fit.
      genfit::MaterialEffects *me = genfit::MaterialEffects::getInstance();
      me->setNoEffects(true);
      ok = doFit();
      me->setNoEffects(false); // restore for the next (clean) track
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
   ft->SetTrackMetadata(std::move(meta));

   return owner.release();
}

} // namespace EventFit
