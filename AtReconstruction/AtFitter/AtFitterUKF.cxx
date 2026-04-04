#include "AtFitterUKF.h"

#include "AtFitMetadata.h"
#include "AtFitTrackMetadata.h"
#include "AtFittedTrack.h"
#include "AtHitCluster.h"
#include "AtKinematics.h"
#include "AtPropagator.h"
#include "AtTrack.h"
#include "AtTrackTransformer.h"

#include <FairLogger.h>

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <TMath.h>
#include <TMatrixD.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <csignal>
#include <csetjmp>
#include <memory>

// Signal handler for catching segfaults during UKF fitting
namespace {
static jmp_buf sJmpBuf;
static bool sHandlerActive = false;
static struct sigaction sOldHandler;

void ukfSegvHandler(int sig)
{
   if (sHandlerActive)
      longjmp(sJmpBuf, 1);
   // If handler not active, call the previous handler
   if (sOldHandler.sa_handler)
      sOldHandler.sa_handler(sig);
}
} // namespace
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

   // --- 1. Get clusters (already ordered by PRA::OrderClustersAlongTrack) ---
   auto *clusters = track->GetHitClusterArray();

   // Convert digi→lab coordinates if ZPadPlane is set.
   // PRA orders clusters along the arc starting from the vertex (highest Z_digi).
   // After Z flip the arc order is preserved but we reverse so that vertex
   // (now lowest Z_lab) comes first — matching the UKF propagation direction.
   if (fZPadPlane > 0) {
      for (auto &cl : *clusters)
         cl.SetPosition({cl.GetPosition().X(), cl.GetPosition().Y(), fZPadPlane - cl.GetPosition().Z()});
   }

   // Filter: minimum spacing + trim Bragg peak end
   {
      std::vector<AtHitCluster> filtered;
      if (!clusters->empty())
         filtered.push_back(clusters->front());
      for (size_t i = 1; i < clusters->size(); i++) {
         double dist = (clusters->at(i).GetPosition() - filtered.back().GetPosition()).R();
         if (dist >= fMinClusterSpacing)
            filtered.push_back(clusters->at(i));
      }
      // Trim last 10% (near Bragg peak where proton stops)
      int nKeep = std::max(fMinClusters, static_cast<int>(filtered.size() * 0.9));
      if (static_cast<int>(filtered.size()) > nKeep)
         filtered.resize(nKeep);
      *clusters = std::move(filtered);
   }

   if (static_cast<int>(clusters->size()) < fMinClusters) {
      LOG(debug) << "AtFitterUKF: " << clusters->size() << " clusters < " << fMinClusters << ". Skipping.";
      return nullptr;
   }

   // Skip tracks with invalid geometry (NaN from failed RANSAC fit)
   if (std::isnan(track->GetGeoTheta()) || std::isnan(track->GetGeoRadius()) || track->GetGeoRadius() <= 0) {
      LOG(debug) << "AtFitterUKF: track has invalid geometry (NaN theta or radius). Skipping.";
      return nullptr;
   }

   // Skip beam-like tracks (lab theta < fMinLabTheta)
   // GeoTheta is in digi frame: theta_lab = 180 - theta_digi
   double thetaLab = 180.0 - track->GetGeoTheta() * 180.0 / M_PI;
   if (thetaLab < fMinLabTheta || thetaLab > (180.0 - fMinLabTheta)) {
      LOG(debug) << "AtFitterUKF: track theta_lab=" << thetaLab << " deg, below threshold " << fMinLabTheta
                 << ". Skipping (beam-like).";
      return nullptr;
   }

   // --- 1b. Adaptive re-clustering based on Brho momentum estimate ---
   // Low-energy tracks benefit from larger cluster radius (more averaging).
   // Re-cluster BEFORE Z conversion (needs digi-frame coordinates).
   if (fAdaptiveClustering) {
      double pEst = GetInitialMomentum(track).R(); // Brho estimate
      if (fMomentumSeed > 0)
         pEst = fMomentumSeed;
      double keEst = std::sqrt(pEst * pEst + fMass_MeV * fMass_MeV) - fMass_MeV;

      // Choose radius/distance based on estimated KE
      double radius, distance;
      if (keEst < 3.0) {
         radius = 25.0;
         distance = 15.0; // Large radius for low-energy (short) tracks
      } else if (keEst < 8.0) {
         radius = 20.0;
         distance = 15.0; // Standard overlapping
      } else {
         radius = 20.0;
         distance = 15.0; // Same for high energy — already good
      }

      // Re-cluster with adapted parameters (works in digi frame)
      // Need to undo Z conversion if already applied — but we haven't yet
      track->ResetHitClusterArray();
      AtTools::AtTrackTransformer transformer;
      transformer.ClusterizeSmooth3D(*track, radius, distance);

      // Re-order clusters along the track
      // (simplified: just use the existing cluster order from ClusterizeSmooth3D)
      // The PRA's OrderClustersAlongTrack would need the track to be in the PRA context
      // so we do a simple highest-Z-first ordering here
      auto *newCl = track->GetHitClusterArray();
      if (newCl->size() >= 2) {
         int seedIdx = 0;
         for (int ci = 1; ci < (int)newCl->size(); ci++) {
            if (newCl->at(ci).GetPosition().Z() > newCl->at(seedIdx).GetPosition().Z())
               seedIdx = ci;
         }
         // NN walk from seed
         std::vector<int> order;
         std::vector<bool> used(newCl->size(), false);
         order.push_back(seedIdx);
         used[seedIdx] = true;
         for (int step = 1; step < (int)newCl->size(); step++) {
            auto cur = newCl->at(order.back()).GetPosition();
            double bestD = 1e9;
            int bestI = -1;
            for (int ci = 0; ci < (int)newCl->size(); ci++) {
               if (used[ci])
                  continue;
               double dd = (newCl->at(ci).GetPosition() - cur).R();
               if (dd < bestD) {
                  bestD = dd;
                  bestI = ci;
               }
            }
            if (bestI < 0)
               break;
            order.push_back(bestI);
            used[bestI] = true;
         }
         std::vector<AtHitCluster> ordered;
         for (int idx : order)
            ordered.push_back(newCl->at(idx));
         *newCl = std::move(ordered);
      }

      clusters = track->GetHitClusterArray(); // refresh pointer
      LOG(debug) << "Adaptive clustering: KE_est=" << keEst << " MeV, r=" << radius << " d=" << distance
                 << " → " << clusters->size() << " clusters";
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

   // --- Iterative fitting ---
   // On iteration 0: use Brho seed with fMomSigmaFrac covariance
   // On subsequent iterations: use previous result as seed with tighter covariance
   double momSigmaFrac = fMomSigmaFrac;
   bool fitConverged = true;
   for (int iter = 0; iter < fNIterations; iter++) {
      fitConverged = true;
      if (iter > 0) {
         // Use the smoothed vertex state from the previous iteration as seed
         auto &smoothedStates = fUKF->GetSmoothedStates();
         if (!smoothedStates.empty()) {
            auto &s0 = smoothedStates[0];
            double pPrev = s0[3];
            double thPrev = s0[4];
            double phPrev = s0[5];
            if (pPrev > 0 && !std::isnan(pPrev)) {
               ROOT::Math::Polar3DVector momPrev(pPrev, thPrev, phPrev);
               initialMom = ROOT::Math::XYZVector(momPrev);
               initialPos = clusters->front().GetPosition();
               p_MeV = pPrev;
               momSigmaFrac = fMomSigmaFrac; // Same sigma on subsequent iterations
               LOG(debug) << "AtFitterUKF: iter " << iter << " seed p=" << pPrev << " MeV/c";
            }
         }
      }

   LOG(debug) << "AtFitterUKF: iter " << iter << "/" << fNIterations << " seed p=" << p_MeV << " MeV/c";

   // --- 3. Set initial state (also calls Reset internally) ---
   {
      // Use current momSigmaFrac for this iteration's covariance
      double savedFrac = fMomSigmaFrac;
      fMomSigmaFrac = (iter == 0) ? fMomSigmaFrac : momSigmaFrac;
      fUKF->SetInitialState(initialPos, initialMom, GetInitialCovariance(p_MeV));
      fMomSigmaFrac = savedFrac;
   }
   fUKF->SetMeasCov(GetMeasCovariance());

   // --- 4. Forward filter pass ---
   auto fitStart = std::chrono::steady_clock::now();

   // Install signal handler to catch segfaults from bad sigma points
   struct sigaction sa;
   sa.sa_handler = ukfSegvHandler;
   sigemptyset(&sa.sa_mask);
   sa.sa_flags = 0;
   sigaction(SIGSEGV, &sa, &sOldHandler);

   if (setjmp(sJmpBuf) != 0) {
      // We get here if a segfault was caught
      LOG(warn) << "AtFitterUKF: segfault caught during fit of track " << track->GetTrackID() << ", skipping";
      sHandlerActive = false;
      sigaction(SIGSEGV, &sOldHandler, nullptr); // Restore original handler
      return nullptr;
   }
   sHandlerActive = true;

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

         // Validate predicted state: check for NaN, Inf, and extreme values
         auto &state = fUKF->vecX();
         auto &cov = fUKF->matP();
         bool stateValid = true;
         for (int d = 0; d < 6; d++) {
            if (std::isnan(state[d]) || std::isinf(state[d]) || std::abs(state[d]) > 1e6) {
               stateValid = false;
               break;
            }
         }
         for (int d = 0; d < 6 && stateValid; d++) {
            if (std::isnan(cov(d, d)) || std::isinf(cov(d, d)) || cov(d, d) < 0 || cov(d, d) > 1e10) {
               stateValid = false;
            }
         }
         if (!stateValid) {
            LOG(debug) << "AtFitterUKF: invalid state/covariance at cluster " << i << ", aborting fit";
            fitConverged = false;
            break;
         }

         fUKF->correctUKF(meas);
      }
   } catch (const std::exception &e) {
      LOG(warn) << "AtFitterUKF: forward pass failed for track " << track->GetTrackID() << ": " << e.what();
      fitConverged = false;
   }

   // Disable signal handler
   sHandlerActive = false;
   sigaction(SIGSEGV, &sOldHandler, nullptr);

   // --- 5. RTS smoother ---
   if (fitConverged) {
      try {
         fUKF->smoothUKF();
      } catch (const std::exception &e) {
         LOG(warn) << "AtFitterUKF: smoother failed for track " << track->GetTrackID() << ": " << e.what();
         fitConverged = false;
      }
   }

   if (!fitConverged)
      break; // Don't iterate if this pass failed

   } // End of iteration loop

   // --- 6. Extract results ---
   const auto &smoothedStates = fUKF->GetSmoothedStates();

   // Build vertex kinematics from the first smoothed state.
   // smoothedStates[0] is at the first cluster, not the true vertex.
   double vx = 0, vy = 0, vz = 0, p_s = p_MeV, theta_s = 0, phi_s = 0;
   if (fitConverged && !smoothedStates.empty()) {
      const auto &s0 = smoothedStates[0];
      vx = s0[0];
      vy = s0[1];
      vz = s0[2];
      p_s = s0[3];
      theta_s = s0[4];
      phi_s = s0[5];

      // Convert angles from Z-flipped frame back to lab frame
      if (fZPadPlane > 0) {
         theta_s = M_PI - theta_s; // Z flip reverses theta
         phi_s = phi_s + M_PI;     // Flip phi by 180°
         if (phi_s > M_PI)
            phi_s -= 2 * M_PI; // Keep in [-π, π]
      }

      // --- Back-extrapolation to beam axis ---
      // The first cluster is some distance from the true vertex at (0,0).
      // Use the actual first cluster position (not the smoothed state which
      // may have been pulled by the smoother).
      auto firstClPos = clusters->front().GetPosition();
      double rXY = std::sqrt(firstClPos.X() * firstClPos.X() + firstClPos.Y() * firstClPos.Y());
      if (rXY > 1.0 && p_s > 0) { // Only extrapolate if meaningfully off-axis
         // Approximate path length from vertex to first cluster:
         // The track curves, so the path is longer than the straight-line distance.
         // Use rXY / sin(theta) as an estimate of the arc length.
         double sinTheta = std::sin(theta_s);
         double pathLength = (sinTheta > 0.1) ? rXY / sinTheta : rXY;

         // Get energy loss per mm at the current momentum
         double KE_at_cluster = std::sqrt(p_s * p_s + fMass_MeV * fMass_MeV) - fMass_MeV;
         if (auto *elossModel = fUKF->GetPropagator().GetELossModel()) {
            double dEdx = elossModel->GetdEdx(KE_at_cluster); // MeV/mm
            double eLost = dEdx * pathLength;                  // Energy lost between vertex and first cluster

            // The proton had MORE energy at the vertex
            double KE_at_vertex = KE_at_cluster + eLost;
            double p_at_vertex = std::sqrt(KE_at_vertex * KE_at_vertex + 2 * KE_at_vertex * fMass_MeV);

            LOG(debug) << "Back-extrapolation: rXY=" << rXY << "mm, path=" << pathLength
                       << "mm, dEdx=" << dEdx << ", eLost=" << eLost
                       << "MeV, p: " << p_s << "→" << p_at_vertex << " MeV/c";

            p_s = p_at_vertex;
         }

         // Extrapolate position back to beam axis along the momentum direction
         ROOT::Math::Polar3DVector momDir(1.0, theta_s, phi_s);
         ROOT::Math::XYZVector dir(momDir);
         // Project back: how far along -dir to reach x²+y²≈0
         // Approximate: move along -dir by pathLength
         vx -= dir.X() * pathLength;
         vy -= dir.Y() * pathLength;
         vz -= dir.Z() * pathLength;
      }
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
