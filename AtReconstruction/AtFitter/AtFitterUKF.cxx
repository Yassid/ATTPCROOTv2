#include "AtFitterUKF.h"

#include "AtFitMetadata.h"
#include "AtFitTrackMetadata.h"
#include "AtFittedTrack.h"
#include "AtHitCluster.h"
#include "AtKinematics.h"
#include "AtPropagator.h"
#include "AtTrack.h"

#include <FairLogger.h>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <TMath.h>
#include <TMatrixD.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include "kalman_filter/TrackFitterUKF.h"

namespace EventFit {

// Destructor defined here so that kf::TrackFitterUKF is a complete type
// (its header is included above via TrackFitterUKF.h).
AtFitterUKF::~AtFitterUKF() = default;

namespace {
// Elementary charge in Coulombs — used to derive integer charge number Z.
constexpr double kElectronCharge = 1.602176634e-19;
// Proton mass in amu equivalent — used for mass unit conversion.
constexpr double kAMU_MeV = 931.494;
// Initial angular uncertainty in the state covariance (degrees, converted to radians when squared).
constexpr double kInitAngSigma_deg = 1.0;
} // namespace

AtFitterUKF::AtFitterUKF(double charge, double mass_MeV, std::unique_ptr<AtTools::AtELossModel> elossModel)
   : fCharge(charge), fMass_MeV(mass_MeV), fELossModel(std::move(elossModel))
{
}

void AtFitterUKF::SetBField(ROOT::Math::XYZVector bField)
{
   fBField = bField;
   if (fUKF)
      fUKF->SetBField(bField);
}

void AtFitterUKF::SetEField(ROOT::Math::XYZVector eField)
{
   fEField = eField;
   if (fUKF)
      fUKF->SetEField(eField);
}

void AtFitterUKF::SetUKFParameters(double alpha, double beta, double kappa)
{
   fAlpha = alpha;
   fBeta = beta;
   fKappa = kappa;
   if (fUKF)
      fUKF->setParameters(static_cast<float>(alpha), static_cast<float>(beta), static_cast<float>(kappa));
}

void AtFitterUKF::InitUKF()
{
   AtTools::AtPropagator propagator(fCharge, fMass_MeV, std::move(fELossModel));
   propagator.SetBField(fBField);
   propagator.SetEField(fEField);

   auto stepper = std::make_unique<AtTools::AtRK4Stepper>();
   fUKF = std::make_unique<kf::TrackFitterUKF>(std::move(propagator), std::move(stepper));
   fUKF->setParameters(static_cast<float>(fAlpha), static_cast<float>(fBeta), static_cast<float>(fKappa));
   fUKF->fEnableEnStraggling = fEnableEnStraggling;
}

double AtFitterUKF::GetBrho(AtTrack *track) const
{
   double radius = track->GetGeoRadius() / 1000.0; // mm → m
   double theta = track->GetGeoTheta();
   double sinTheta = std::sin(theta);
   if (std::abs(sinTheta) < 1e-9)
      sinTheta = 1e-9;
   return fBField.Z() * radius / sinTheta; // T·m
}

ROOT::Math::XYZPoint AtFitterUKF::GetInitialPosition(AtTrack *track) const
{
   return track->GetHitClusterArray()->front().GetPosition(); // mm
}

ROOT::Math::XYZVector AtFitterUKF::GetInitialMomentum(AtTrack *track) const
{
   double brho = GetBrho(track);
   int Z = static_cast<int>(std::round(fCharge / kElectronCharge));
   double mass_amu = fMass_MeV / kAMU_MeV;

   auto [p_MeV, ke_MeV] = AtTools::Kinematics::GetMomFromBrho(mass_amu, Z, brho);

   double theta = track->GetGeoTheta();
   double phi = track->GetGeoPhi();
   ROOT::Math::Polar3DVector momPolar(p_MeV, theta, phi);
   return {momPolar.X(), momPolar.Y(), momPolar.Z()};
}

TMatrixD AtFitterUKF::GetInitialCovariance(double p_mag_MeV) const
{
   TMatrixD cov(6, 6);
   cov.Zero();
   double sigma_pos2 = fMeasSigma_mm * fMeasSigma_mm;
   double sigma_mom2 = (fMomSigmaFrac * p_mag_MeV) * (fMomSigmaFrac * p_mag_MeV);
   double sigma_ang2 = (kInitAngSigma_deg * TMath::Pi() / 180.0) * (kInitAngSigma_deg * TMath::Pi() / 180.0);
   cov(0, 0) = sigma_pos2;
   cov(1, 1) = sigma_pos2;
   cov(2, 2) = sigma_pos2;
   cov(3, 3) = sigma_mom2;
   cov(4, 4) = sigma_ang2;
   cov(5, 5) = sigma_ang2;
   return cov;
}

TMatrixD AtFitterUKF::GetMeasCovariance() const
{
   TMatrixD cov(3, 3);
   cov.Zero();
   double sigma2 = fMeasSigma_mm * fMeasSigma_mm;
   cov(0, 0) = sigma2;
   cov(1, 1) = sigma2;
   cov(2, 2) = sigma2;
   return cov;
}

AtFittedTrack *
AtFitterUKF::GetFittedTrack(AtTrack *track, AtFitMetadata *fitMetadata, AtRawEvent * /*rawEvent*/, AtEvent * /*event*/)
{
   // Lazy initialisation of the UKF on first call.
   if (!fUKF)
      InitUKF();

   // --- 1. Get clusters (PRA already orders them by Z/time along the track) ---
   auto *clusters = track->GetHitClusterArray();

   // Convert digi→lab coordinates if ZPadPlane is set
   if (fZPadPlane > 0) {
      for (auto &cl : *clusters) {
         auto pos = cl.GetPosition();
         cl.SetPosition({pos.X(), pos.Y(), fZPadPlane - pos.Z()});
      }
      // Reverse: PRA ordered by ascending Z_digi (Bragg peak → vertex),
      // after flip this becomes descending Z_lab. Reverse to get vertex first.
      std::reverse(clusters->begin(), clusters->end());
   }
   if (static_cast<int>(clusters->size()) < fMinClusters) {
      LOG(info) << "AtFitterUKF: track " << track->GetTrackID() << " has " << clusters->size()
                << " clusters, fewer than fMinClusters=" << fMinClusters << ". Skipping.";
      return nullptr;
   }

   // --- 2. Momentum seed ---
   ROOT::Math::XYZPoint initialPos = GetInitialPosition(track);
   ROOT::Math::XYZVector initialMom = GetInitialMomentum(track);
   double p_MeV = initialMom.R();

   // Override momentum magnitude if seed is set
   if (fMomentumSeed > 0)
      p_MeV = fMomentumSeed;

   // If we converted to lab frame, the GeoTheta/Phi from pattern recognition
   // are in the digi frame and no longer match the cluster positions.
   // Re-derive the momentum direction from the first two clusters instead.
   if (fZPadPlane > 0 && clusters->size() >= 2) {
      auto dir = clusters->at(1).GetPosition() - clusters->at(0).GetPosition();
      initialMom = p_MeV * dir.Unit();
   }

   LOG(debug) << "AtFitterUKF: seed p=" << p_MeV << " MeV/c, pos=(" << initialPos.X() << "," << initialPos.Y() << ","
              << initialPos.Z() << "), mom=(" << initialMom.X() << "," << initialMom.Y() << "," << initialMom.Z()
              << ")";

   // --- 3. Set initial state (also calls Reset internally) ---
   fUKF->SetInitialState(initialPos, initialMom, GetInitialCovariance(p_MeV));
   fUKF->SetMeasCov(GetMeasCovariance());

   // --- 4. Forward filter pass ---
   bool fitConverged = true;
   auto fitStart = std::chrono::steady_clock::now();
   try {
      for (size_t i = 1; i < clusters->size(); ++i) {
         // Timeout: abort if fit takes too long (default 2 seconds per track)
         auto elapsed = std::chrono::steady_clock::now() - fitStart;
         if (std::chrono::duration_cast<std::chrono::milliseconds>(elapsed).count() > fMaxFitTime_ms) {
            LOG(warn) << "AtFitterUKF: timeout after " << fMaxFitTime_ms << " ms at cluster " << i << "/"
                      << clusters->size();
            fitConverged = false;
            break;
         }

         ROOT::Math::XYZPoint meas = clusters->at(i).GetPosition();

         // Update measurement covariance from cluster if enabled
         if (fUsePerClusterCov) {
            const auto &clCov = clusters->at(i).GetCovMatrix();
            TMatrixD measCov(3, 3);
            for (int r = 0; r < 3; ++r)
               for (int c = 0; c < 3; ++c)
                  measCov(r, c) = clCov(r, c);
            double minVar = 0.01;
            for (int d = 0; d < 3; ++d)
               if (measCov(d, d) < minVar)
                  measCov(d, d) = minVar;
            fUKF->SetMeasCov(measCov);
         }

         fUKF->predictUKF(meas);
         fUKF->correctUKF(meas);
      }
   } catch (const std::exception &e) {
      LOG(warn) << "AtFitterUKF: forward pass failed for track " << track->GetTrackID() << ": " << e.what();
      fitConverged = false;
   }

   // --- 5. RTS smoother ---
   if (fitConverged) {
      try {
         fUKF->smoothUKF();
      } catch (const std::exception &e) {
         LOG(warn) << "AtFitterUKF: smoother failed for track " << track->GetTrackID() << ": " << e.what();
         fitConverged = false;
      }
   }

   // --- 6. Extract results ---
   const auto &smoothedStates = fUKF->GetSmoothedStates();

   // Build vertex kinematics from the first smoothed state.
   // smoothedStates[0] is the vertex: SetInitialState seeds index 0, and the RTS smoother
   // back-propagates all the way to index 0, so this is the best-estimate vertex state.
   double vx = 0, vy = 0, vz = 0, p_s = p_MeV, theta_s = 0, phi_s = 0;
   if (fitConverged && !smoothedStates.empty()) {
      const auto &s0 = smoothedStates[0];
      vx = s0[0];
      vy = s0[1];
      vz = s0[2];
      p_s = s0[3];
      theta_s = s0[4];
      phi_s = s0[5];
   }

   double KE = std::sqrt(p_s * p_s + fMass_MeV * fMass_MeV) - fMass_MeV;

   // --- 7. Build AtFittedTrack ---
   // Use a unique_ptr internally so that any exception thrown by the setters below does not leak.
   auto fittedTrackOwner = std::make_unique<AtFittedTrack>();
   auto *fittedTrack = fittedTrackOwner.get();
   fittedTrack->SetTrackID(track->GetTrackID());
   fittedTrack->SetKinematics(KE, theta_s, phi_s);
   fittedTrack->SetVertex(ROOT::Math::XYZVector(vx, vy, vz));

   int Z = static_cast<int>(std::round(fCharge / kElectronCharge));
   double mass_amu = fMass_MeV / kAMU_MeV;
   fittedTrack->SetParticleInfo(std::to_string(Z), Z, mass_amu);

   // Smoothed positions: indices 1..N-1 (skip the vertex at index 0).
   std::vector<ROOT::Math::XYZPoint> smoothedPos;
   if (fitConverged && smoothedStates.size() > 1) {
      smoothedPos.reserve(smoothedStates.size() - 1);
      for (size_t i = 1; i < smoothedStates.size(); ++i) {
         const auto &s = smoothedStates[i];
         smoothedPos.emplace_back(s[0], s[1], s[2]);
      }
   }
   fittedTrack->SetSmoothedPositions(std::move(smoothedPos));

   // Track length from consecutive smoothed positions.
   double trackLength = 0;
   const auto &spos = fittedTrack->GetSmoothedPositions();
   if (!spos.empty()) {
      ROOT::Math::XYZPoint prev(vx, vy, vz);
      for (const auto &pt : spos) {
         trackLength += (pt - prev).R();
         prev = pt;
      }
   }

   // Chi2: sum of squared distances between smoothed positions and measurements / ndf.
   double chi2 = 0;
   int nClusters = static_cast<int>(clusters->size());
   if (fitConverged && static_cast<int>(smoothedStates.size()) == nClusters) {
      for (int i = 1; i < nClusters; ++i) {
         const auto &s = smoothedStates[i];
         ROOT::Math::XYZPoint smoothedPt(s[0], s[1], s[2]);
         ROOT::Math::XYZPoint measPt = clusters->at(i).GetPosition();
         double dist = (smoothedPt - measPt).R();
         chi2 += dist * dist / (fMeasSigma_mm * fMeasSigma_mm);
      }
   }
   // Subtract the number of fitted parameters (6: x,y,z,p,theta,phi) from the
   // degrees of freedom.  Each cluster contributes 3 measurements (x,y,z), so
   // the total number of measurements is 3*(nClusters-1) (we skip the seed point).
   // ndf = nMeasurements - nParams = 3*(nClusters-1) - 6.
   int ndf = std::max(1, 3 * (nClusters - 1) - 6);

   // Per-track metadata.
   auto trackMeta = std::make_unique<AtFitTrackMetadata>();
   trackMeta->SetTrackID(track->GetTrackID());
   trackMeta->SetFitConverged(fitConverged);
   trackMeta->SetChi2(chi2 / ndf);
   trackMeta->SetNdf(ndf);
   fittedTrack->SetTrackMetadata(std::move(trackMeta));

   // Track properties.
   fittedTrack->SetTrackPropertiesStruct(ROOT::Math::XYZVector(initialPos.X(), initialPos.Y(), initialPos.Z()),
                                         ROOT::Math::XYZVector(vx, vy, vz), /*extrapolatedDistance=*/0,
                                         /*distancePOCA=*/0, trackLength, trackLength, /*estimateTotalCharge=*/0,
                                         /*trackPoints=*/nClusters);

   // Optionally record metadata in the event-level AtFitMetadata.
   if (fitMetadata) {
      AtFitMetadata::TrackMetadatasVector metaVec;
      auto meta2 = std::make_unique<AtFitTrackMetadata>();
      meta2->SetTrackID(track->GetTrackID());
      meta2->SetFitConverged(fitConverged);
      meta2->SetChi2(chi2 / ndf);
      meta2->SetNdf(ndf);
      metaVec.push_back(std::move(meta2));
      fitMetadata->SetTrackMetadatasVector(track->GetTrackID(), std::move(metaVec));
   }

   return fittedTrackOwner.release();
}

} // namespace EventFit
