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
   fUKF->fELossScaleFactor = fELossScaleFactor;
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

   // Skip tracks with invalid geometry (NaN from failed RANSAC fit)
   {
      const bool badRadius = std::isnan(track->GetGeoRadius()) || track->GetGeoRadius() <= 0;
      const bool nanTheta = std::isnan(track->GetGeoTheta());
      if (badRadius) {
         LOG(info) << "AtFitterUKF DROP[bad-geom-radius]: R=" << track->GetGeoRadius() << ". Skipping.";
         return nullptr;
      }
      // PRA sometimes returns NaN GeoTheta with a valid radius (degenerate
      // (arc, z) regression). Recover theta from the cluster sequence so the
      // track isn't discarded.
      if (nanTheta) {
         // PRA sometimes returns NaN GeoTheta when the (arc, z) line fit is
         // degenerate. Recover θ from the track's endpoint geometry. Prefer
         // the cluster array if populated; otherwise fall back to the raw
         // hit array (always present after PRA). This avoids dropping tracks
         // before adaptive re-clustering has had a chance to populate the
         // cluster array further down the flow.
         auto *cls0 = track->GetHitClusterArray();
         double x0 = 0, y0 = 0, z0 = 0, xN = 0, yN = 0, zN = 0;
         bool haveEndpoints = false;
         if (cls0 != nullptr && cls0->size() >= 2) {
            const auto &p0 = cls0->front().GetPosition();
            const auto &pN = cls0->back().GetPosition();
            x0 = p0.X(); y0 = p0.Y(); z0 = p0.Z();
            xN = pN.X(); yN = pN.Y(); zN = pN.Z();
            haveEndpoints = true;
         } else {
            const auto &hits = track->GetHitArray();
            if (hits.size() >= 2) {
               const auto &p0 = hits.front()->GetPosition();
               const auto &pN = hits.back()->GetPosition();
               x0 = p0.X(); y0 = p0.Y(); z0 = p0.Z();
               xN = pN.X(); yN = pN.Y(); zN = pN.Z();
               haveEndpoints = true;
            }
         }
         if (!haveEndpoints) {
            LOG(info) << "AtFitterUKF DROP[bad-geom-noendpoints]: NaN theta and "
                         "<2 clusters/hits. Skipping.";
            return nullptr;
         }
         const double dx = xN - x0;
         const double dy = yN - y0;
         const double dz = zN - z0;
         const double r3 = std::sqrt(dx * dx + dy * dy + dz * dz);
         if (r3 < 1e-3) {
            LOG(info) << "AtFitterUKF DROP[bad-geom-coincide]: NaN theta + endpoints coincide. Skipping.";
            return nullptr;
         }
         const double thetaFallback = std::acos(dz / r3);
         track->SetGeoTheta(thetaFallback);
         if (std::isnan(track->GetGeoPhi()))
            track->SetGeoPhi(std::atan2(dy, dx));
         LOG(info) << "AtFitterUKF: PRA theta NaN, recovered: "
                   << thetaFallback * 180. / M_PI << " deg (from "
                   << (cls0 && cls0->size() >= 2 ? "clusters" : "hits") << ")";
      }
   }

   // Skip beam-like tracks (lab theta < fMinLabTheta)
   // GeoTheta is in digi frame: theta_lab = 180 - theta_digi
   double thetaLab = 180.0 - track->GetGeoTheta() * 180.0 / M_PI;
   if (thetaLab < fMinLabTheta || thetaLab > (180.0 - fMinLabTheta)) {
      LOG(info) << "AtFitterUKF DROP[beam-like]: theta_lab=" << thetaLab << " deg, threshold " << fMinLabTheta
                << ". Skipping.";
      return nullptr;
   }

   // --- 1b. Adaptive re-clustering based on Brho momentum estimate ---
   // Low-energy tracks benefit from larger cluster radius (more averaging).
   // Re-clustering works in digi frame (hits are never Z-flipped).
   if (fAdaptiveClustering) {
      double pEst = GetInitialMomentum(track).R(); // Brho estimate
      if (fMomentumSeed > 0)
         pEst = fMomentumSeed;
      double keEst = std::sqrt(pEst * pEst + fMass_MeV * fMass_MeV) - fMass_MeV;

      // Choose radius/distance:
      //   1. Hard override via SetClusteringParams (highest priority)
      //   2. Arc-length-aware mode via SetTargetClusters (PUMA-style)
      //   3. Legacy KE-keyed defaults
      double radius, distance;
      if (fClusterRadiusOverride > 0 && fClusterDistanceOverride > 0) {
         radius = fClusterRadiusOverride;
         distance = fClusterDistanceOverride;
      } else if (fTargetClusters > 0) {
         // Compute arc length from the PRA circle: R · |Δφ| spanned by hits
         // around GeoCenter. Falls back to 3D extent of the hit array if the
         // circle parameters look bad.
         const auto &hits = track->GetHitArray();
         double arc = 0;
         auto cen = track->GetGeoCenter();
         double cx = cen.first;
         double cy = cen.second;
         double R = track->GetGeoRadius();
         bool circleValid = (hits.size() >= 2) && (R > 1.0) && std::isfinite(R);
         if (circleValid) {
            double phiMin = std::numeric_limits<double>::infinity();
            double phiMax = -std::numeric_limits<double>::infinity();
            for (const auto &h : hits) {
               const auto &p = h->GetPosition();
               double phi = std::atan2(p.Y() - cy, p.X() - cx);
               if (phi < phiMin) phiMin = phi;
               if (phi > phiMax) phiMax = phi;
            }
            // Span doesn't handle wrap-around but the PRA circle for PUMA
            // arcs that don't loop is reliable here.
            double dphi = phiMax - phiMin;
            arc = R * dphi;
         }
         if (!(arc > 1.0)) {
            // Fall back to first-last 3D distance — not ideal for curved
            // tracks but better than nothing.
            if (hits.size() >= 2) {
               const auto &p0 = hits.front()->GetPosition();
               const auto &pN = hits.back()->GetPosition();
               double dx = pN.X() - p0.X();
               double dy = pN.Y() - p0.Y();
               double dz = pN.Z() - p0.Z();
               arc = std::sqrt(dx * dx + dy * dy + dz * dz);
            }
         }
         distance = std::max(fAdaptiveDistMin,
                             std::min(fAdaptiveDistMax, arc / fTargetClusters));
         radius = 1.5 * distance;
         LOG(debug) << "Arc-length adaptive clustering: arc=" << arc
                    << " mm, target=" << fTargetClusters
                    << " → radius=" << radius << " distance=" << distance << " mm";
      } else if (keEst < 3.0) {
         radius = 25.0;
         distance = 15.0; // Large radius for low-energy (short) tracks
      } else if (keEst < 8.0) {
         radius = 20.0;
         distance = 15.0; // Standard overlapping
      } else {
         radius = 20.0;
         distance = 15.0; // Same for high energy — already good
      }

      // Re-cluster with adapted parameters (works in digi frame — hits are never Z-flipped).
      // The digi→lab Z conversion is applied after this block.
      track->ResetHitClusterArray();
      AtTools::AtTrackTransformer transformer;
      if (fUseArcWalk) {
         // Geometry-based arc ordering, gap-immune. Target = fTargetClusters
         // when set (PUMA path); otherwise fall back to the ArcWalk default.
         const int target = (fTargetClusters > 0) ? fTargetClusters : 25;
         transformer.ClusterizeArcWalk(*track, target, fArcWalkMinHits, fArcWalkKNN);
      } else {
         transformer.ClusterizeSmooth3D(*track, radius, distance);
      }

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

   // Filter: minimum spacing + trim Bragg peak end.
   // Applied after adaptive re-clustering so it operates on the final cluster set.
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
      LOG(info) << "AtFitterUKF DROP[few-clusters]: " << clusters->size() << " < " << fMinClusters << ". Skipping.";
      return nullptr;
   }

   // Convert digi→lab coordinates if ZPadPlane is set.
   // This must happen AFTER adaptive re-clustering (which produces digi-frame clusters).
   // PRA orders clusters along the arc starting from the vertex (highest Z_digi).
   // After Z flip the arc order is preserved; vertex (now lowest Z_lab) comes first.
   if (fZPadPlane > 0) {
      for (auto &cl : *clusters)
         cl.SetPosition({cl.GetPosition().X(), cl.GetPosition().Y(), fZPadPlane - cl.GetPosition().Z()});
   }

   // --- 2. Momentum seed ---
   ROOT::Math::XYZPoint initialPos = GetInitialPosition(track);
   ROOT::Math::XYZVector initialMom = GetInitialMomentum(track);
   double p_MeV = initialMom.R();

   // Override momentum magnitude if seed is set
   if (fMomentumSeed > 0)
      p_MeV = fMomentumSeed;

   // Override the GeoTheta/GeoPhi-derived direction when:
   //   * fZPadPlane > 0 — lab-frame conversion makes the digi-frame Geo
   //     angles inconsistent with the cluster positions (legacy 16C+p path),
   //   * fUseClusterDirSeed — PRA's (arclength,z) line fit gives an
   //     ambiguous-half-sphere theta (e.g. PUMA upward tracks see
   //     GeoTheta=115° instead of 65°).
   //
   // Two-step seed for fUseClusterDirSeed (works for both PUMA's centred
   // vertex and 16C+p's high-Z vertex):
   //   1. Re-sort clusters by xy-distance from the beam axis (closest first)
   //      so cluster[0] always sits at the vertex end of the track.
   //   2. Initial direction = tangent to the PRA circle at cluster[0],
   //      i.e. perpendicular to (cluster[0] − GeoCenter) in (x,y); the
   //      tangent's sign is disambiguated by aligning with the cluster
   //      chord (no reliance on PRA chargeSign whose B-frame convention
   //      can be wrong). The z-slope of the chord sets the helix pitch.
   //   This avoids the chord-vs-tangent bias on weakly-curved tracks
   //   (e.g. high-momentum PUMA K+) that flipped K-vs-π mass ID.
   if ((fZPadPlane > 0 || fUseClusterDirSeed) && clusters->size() >= 2) {
      if (fUseClusterDirSeed) {
         std::sort(clusters->begin(), clusters->end(),
                   [](const AtHitCluster &a, const AtHitCluster &b) {
                      const auto &pa = a.GetPosition();
                      const auto &pb = b.GetPosition();
                      return (pa.X() * pa.X() + pa.Y() * pa.Y())
                           < (pb.X() * pb.X() + pb.Y() * pb.Y());
                   });
         initialPos = clusters->front().GetPosition();

         const auto &p0 = clusters->at(0).GetPosition();
         const auto &p1 = clusters->at(1).GetPosition();
         const auto cen = track->GetGeoCenter();
         const double R = track->GetGeoRadius();

         double ux = p1.X() - p0.X();
         double uy = p1.Y() - p0.Y();
         double uz = p1.Z() - p0.Z();

         // Use circle tangent for (x,y) direction when PRA gave a sane
         // circle; fall back to the raw chord otherwise.
         if (std::isfinite(R) && R > 0.5) {
            const double rx = p0.X() - cen.first;
            const double ry = p0.Y() - cen.second;
            // Tangent candidate at p0: 90° CCW rotation of the radial.
            double tx = -ry;
            double ty = rx;
            const double tn = std::sqrt(tx * tx + ty * ty);
            if (tn > 0) {
               tx /= tn;
               ty /= tn;
               // Disambiguate sign by aligning with chord direction.
               const double chord_xy_along_t = (p1.X() - p0.X()) * tx
                                             + (p1.Y() - p0.Y()) * ty;
               if (chord_xy_along_t < 0) {
                  tx = -tx;
                  ty = -ty;
               }
               // 3D direction: tangent in xy, chord-derived slope in z.
               const double a = std::abs(chord_xy_along_t);
               const double b = uz;
               const double norm3 = std::sqrt(a * a + b * b);
               if (norm3 > 0) {
                  ux = (a / norm3) * tx;
                  uy = (a / norm3) * ty;
                  uz = b / norm3;
               }
            }
         }

         const double mag = std::sqrt(ux * ux + uy * uy + uz * uz);
         if (mag > 0)
            initialMom = ROOT::Math::XYZVector(p_MeV * ux / mag,
                                               p_MeV * uy / mag,
                                               p_MeV * uz / mag);
      } else {
         // Legacy 16C+p path: raw chord direction.
         auto dir = clusters->at(1).GetPosition() - clusters->at(0).GetPosition();
         initialMom = p_MeV * dir.Unit();
      }
   }

   // --- Iterative fitting ---
   // On iteration 0: use Brho seed with fMomSigmaFrac covariance
   // On subsequent iterations: use previous result as seed with tighter covariance
   double momSigmaFrac = fMomSigmaFrac;
   bool fitConverged = true;
   for (int iter = 0; iter < fNIterations; iter++) {
      fitConverged = true;
      // Reference-track warm start: seed covariance taken from the previous
      // iteration's smoothed vertex (tight) instead of the loose Brho prior.
      TMatrixD warmSeedCov(6, 6);
      bool useWarmSeed = false;
      // Full per-cluster reference-track: previous iteration's smoothed trajectory.
      std::vector<std::array<double, 6>> refTraj;
      bool useRefThisIter = false;
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

               // Warm-start covariance = DIAGONAL of the inflated smoothed-vertex
               // covariance. We deliberately drop the off-diagonal correlations:
               // the full smoothed matrix is near-singular and breaks the
               // sigma-point Cholesky on the re-pass (mass fit drops). The diagonal
               // carries the tight per-component uncertainty (the warm-start benefit)
               // while staying well-conditioned.
               if (fUseRefTrackWarmStart) {
                  auto &smoothedCovs = fUKF->GetSmoothedCovariances();
                  if (smoothedCovs.size() == smoothedStates.size() && !smoothedCovs.empty()) {
                     const auto &c0 = smoothedCovs[0];
                     warmSeedCov.Zero();
                     bool finite = true;
                     double posFloor = fMeasSigma_mm * fMeasSigma_mm;
                     double angFloor = std::pow(0.5 * TMath::Pi() / 180.0, 2); // (0.5 deg)^2
                     for (int d = 0; d < 6 && finite; ++d) {
                        double v = c0(d, d) * fRefTrackInflation;
                        if (!std::isfinite(v) || v < 0) { finite = false; break; }
                        warmSeedCov(d, d) = v;
                     }
                     if (finite) {
                        for (int d = 0; d < 3; ++d)
                           if (warmSeedCov(d, d) < posFloor) warmSeedCov(d, d) = posFloor;
                        for (int d = 4; d < 6; ++d)
                           if (warmSeedCov(d, d) < angFloor) warmSeedCov(d, d) = angFloor;
                        if (warmSeedCov(3, 3) <= 0) warmSeedCov(3, 3) = std::pow(0.05 * pPrev, 2);
                        useWarmSeed = true;
                     }
                  }
               }

               // Full per-cluster reference-track: capture the previous iteration's
               // smoothed trajectory (one state per cluster) before SetInitialState
               // resets the filter. predictUKFRef linearizes around it at every step.
               if (fUseRefTrack && smoothedStates.size() == clusters->size()) {
                  refTraj.resize(smoothedStates.size());
                  bool ok = true;
                  for (size_t s = 0; s < smoothedStates.size() && ok; ++s)
                     for (int d = 0; d < 6; ++d) {
                        double v = smoothedStates[s][d];
                        if (!std::isfinite(v)) { ok = false; break; }
                        refTraj[s][d] = v;
                     }
                  useRefThisIter = ok;
               }
            }
         }
      }

   LOG(debug) << "AtFitterUKF: iter " << iter << "/" << fNIterations << " seed p=" << p_MeV << " MeV/c";

   // --- 3. Set initial state (also calls Reset internally) ---
   {
      // Use current momSigmaFrac for this iteration's covariance
      double savedFrac = fMomSigmaFrac;
      fMomSigmaFrac = (iter == 0) ? fMomSigmaFrac : momSigmaFrac;
      fUKF->SetInitialState(initialPos, initialMom, useWarmSeed ? warmSeedCov : GetInitialCovariance(p_MeV));
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
      LOG(warn) << "AtFitterUKF DROP[segfault]: caught during fit of track " << track->GetTrackID() << ", skipping";
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

         // Reference-track: linearize transport around the previous smoothed state
         // at this cluster (i-1). The previous-iteration smoothed trajectory can
         // contain unphysical points MID-TRACK (negative/huge momentum) — the
         // baseline only ever uses the vertex point so it never noticed, but
         // ref-linearizing around a garbage reference poisons the pass. Validate
         // each reference point and fall back to the standard predict when bad.
         bool refOk = useRefThisIter && (i - 1) < refTraj.size();
         if (refOk) {
            const auto &r = refTraj[i - 1];
            double pp = r[3];
            refOk = (pp > 0 && pp < 1000.0) && std::isfinite(r[0]) && std::isfinite(r[1]) &&
                    std::isfinite(r[2]) && std::isfinite(r[4]) && std::isfinite(r[5]);
         }
         if (refOk)
            fUKF->predictUKFRef(meas, refTraj[i - 1]);
         else
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
   double p_first_cluster_MeV = p_MeV; // pre-back-extrap (smoothed state at first cluster)
   if (fitConverged && !smoothedStates.empty()) {
      const auto &s0 = smoothedStates[0];
      vx = s0[0];
      vy = s0[1];
      vz = s0[2];
      p_s = s0[3];
      theta_s = s0[4];
      phi_s = s0[5];
      p_first_cluster_MeV = p_s; // capture before back-extrap mutates p_s

      // The UKF now operates in lab frame (Z_lab = fZPadPlane - Z_digi), so theta/phi
      // from the smoothed state are already in the correct physics frame.
      // No angle conversion needed.

      // --- Back-extrapolation to beam axis ---
      // Two paths:
      //   (a) fUseHelixBackExtrap = true: closed-form POCA on the PRA circle
      //       (x,y) + helix-pitch z. Use this when tracks curve significantly
      //       in (x,y) so the linear tangent approximation fails (e.g. PUMA
      //       at 4 T).
      //   (b) default: linear step along the initial momentum direction by
      //       pathLength = rXY/sinθ, capped at fBackExtrapMaxPath. Stable for
      //       small-arc tracks (e.g. 16C+p protons) and avoids the Bragg-peak
      //       dE/dx instability that killed the older propagator approach.
      auto firstClPos = clusters->front().GetPosition();
      double rXY = std::sqrt(firstClPos.X() * firstClPos.X() + firstClPos.Y() * firstClPos.Y());
      if (rXY > 1.0 && p_s > 0 && fBackExtrapMaxPath > 0.) {
         double pathLength = 0.;

         // PRA circle parameters (for the helix path)
         auto geoCenter = track->GetGeoCenter();
         double cx = geoCenter.first;
         double cy = geoCenter.second;
         double R = track->GetGeoRadius();
         double dCenter = std::sqrt(cx * cx + cy * cy);
         bool circleValid = fUseHelixBackExtrap && (R > 1.0) && (dCenter > 1.0);

         if (circleValid) {
            // POCA on the PRA circle to (0,0). The two extrema on the circle
            // along the line from origin through center are at
            // (cx,cy)·(1 ± R/d). |d - R| < |d + R|, so the closer one is
            // (cx,cy)·(1 - R/d) for any d > 0.
            double f = 1.0 - R / dCenter;
            double pocaX = cx * f;
            double pocaY = cy * f;

            // Use the raw hit closest to the beam axis as the back-extrap
            // starting point, not the smoothed-state cluster centroid. The
            // centroid sits ~half-cluster-spacing inside the track from the
            // actual first physical hit; on upward-going tracks
            // (vz_first_hit < vz_centroid) this introduces a positive Δz
            // bias of size ~half_cluster · cot(θ). The first raw hit lies
            // on the helix at its true xy/z position, so back-extrap from
            // there along the same arc removes the bias.
            double rxFirst = vx;
            double ryFirst = vy;
            double rzFirst = vz;
            {
               const auto &hits = track->GetHitArray();
               double bestR2 = rxFirst * rxFirst + ryFirst * ryFirst;
               for (const auto &hp : hits) {
                  const auto &p = hp->GetPosition();
                  double r2 = p.X() * p.X() + p.Y() * p.Y();
                  if (r2 < bestR2) {
                     bestR2 = r2;
                     rxFirst = p.X();
                     ryFirst = p.Y();
                     rzFirst = p.Z();
                  }
               }
            }

            // Arc length from raw first hit to POCA along the PRA circle —
            // pick the shorter wrap.
            double phi1 = std::atan2(ryFirst - cy, rxFirst - cx);
            double phi0 = std::atan2(pocaY - cy, pocaX - cx);
            double dPhi = std::abs(phi1 - phi0);
            if (dPhi > M_PI) dPhi = 2 * M_PI - dPhi;
            double arc = R * dPhi;

            // Optional: PRA circle is fit through hits only and typically
            // does NOT pass through the actual beam-axis vertex (offset of
            // a few mm in xy); the POCA-on-circle differs from origin by
            // |dCenter − R|. Extending the arc by the chord from POCA to
            // origin compensates for the corresponding cot(θ) z-bias.
            double endX = pocaX;
            double endY = pocaY;
            if (fForceVertexOnBeamAxis) {
               double chord = std::sqrt(pocaX * pocaX + pocaY * pocaY);
               arc += chord;
               endX = 0;
               endY = 0;
            }
            arc = std::min(arc, fBackExtrapMaxPath);

            // z propagation: helix pitch dz/d(arc_xy) = cot(θ); back-extrap
            // along the helix reverses the sign. θ is constant along the
            // helix so the smoothed-state θ at first cluster is valid here.
            double sinTheta = std::max(std::sin(theta_s), 0.1);
            double cotTheta = std::cos(theta_s) / sinTheta;

            vx = endX;
            vy = endY;
            vz = rzFirst - arc * cotTheta;
            pathLength = arc;

            // Optionally rotate φ_s by the back-extrap arc length / radius so
            // the stored Kinematics.phi is at the vertex, not at the first
            // cluster. Sign: pi- in +Bz curls counterclockwise (dφ/ds > 0),
            // so going BACKWARD in time subtracts |arc/R|; pi+ adds it. θ is
            // unchanged (motion in cyclotron plane is xy, z-component is
            // preserved for uniform Bz).
            if (fUpdateAnglesOnBackExtrap && R > 1.0) {
               const double dphi_mag = arc / R;
               const double signFactor = (fCharge < 0) ? -1.0 : +1.0;
               phi_s += signFactor * dphi_mag;
               // Wrap to (-π, π] to keep downstream consumers happy.
               while (phi_s > M_PI) phi_s -= 2 * M_PI;
               while (phi_s <= -M_PI) phi_s += 2 * M_PI;
            }
         } else {
            // Legacy linear extrapolation
            double sinTheta = std::sin(theta_s);
            pathLength = (sinTheta > 0.1) ? rXY / sinTheta : rXY;
            pathLength = std::min(pathLength, fBackExtrapMaxPath);

            ROOT::Math::Polar3DVector momDir(1.0, theta_s, phi_s);
            ROOT::Math::XYZVector dir(momDir);
            vx -= dir.X() * pathLength;
            vy -= dir.Y() * pathLength;
            vz -= dir.Z() * pathLength;
         }

         // Optional straight-line tail beyond POCA. For setups where the
         // production vertex sits in a field-free region upstream of the
         // chamber, continue from POCA along the back-extrapolated momentum
         // direction (phi_s is at POCA if fUpdateAnglesOnBackExtrap is on)
         // until vx reaches fBackExtrapTargetX. No eloss in the tail
         // (assumed vacuum). Skipped if NaN or already past the target.
         if (!std::isnan(fBackExtrapTargetX) && vx > fBackExtrapTargetX) {
            ROOT::Math::Polar3DVector momDir(1.0, theta_s, phi_s);
            ROOT::Math::XYZVector dir(momDir);
            if (std::abs(dir.X()) > 0.01) {
               double tailPath = (vx - fBackExtrapTargetX) / dir.X();
               vx -= dir.X() * tailPath;
               vy -= dir.Y() * tailPath;
               vz -= dir.Z() * tailPath;
            }
         }

         // Subtract a hardwired digi-z offset if requested. The PSA z formula
         // uses the pulse-peak time bucket without subtracting the
         // peakingTime·vDrift delay that AtPulse added — for PUMA this is a
         // constant +7.5 mm (~+8.6 mm with shape asymmetry) bias on every
         // hit, which propagates straight through the back-extrap to the
         // vertex. Default fVertexZBias_mm = 0 keeps existing experiments
         // unaffected.
         if (fVertexZBias_mm != 0.0)
            vz -= fVertexZBias_mm;

         // Energy correction along the back-extrapolated path. Same in both
         // branches — uses the integrated arc/path length.
         double KE_at_cluster = std::sqrt(p_s * p_s + fMass_MeV * fMass_MeV) - fMass_MeV;
         if (auto *elossModel = fUKF->GetPropagator().GetELossModel()) {
            double dEdx = elossModel->GetdEdx(KE_at_cluster);
            double eLost = dEdx * pathLength;
            double KE_at_vertex = KE_at_cluster + eLost;
            double p_at_vertex = std::sqrt(KE_at_vertex * KE_at_vertex + 2 * KE_at_vertex * fMass_MeV);
            p_s = p_at_vertex;
         }
      }
   }

   double KE = std::sqrt(p_s * p_s + fMass_MeV * fMass_MeV) - fMass_MeV;
   double KE_first_cluster = std::sqrt(p_first_cluster_MeV * p_first_cluster_MeV
                                       + fMass_MeV * fMass_MeV) - fMass_MeV;

   // --- 7. Build AtFittedTrack ---
   // Use a unique_ptr internally so that any exception thrown by the setters below does not leak.
   auto fittedTrackOwner = std::make_unique<AtFittedTrack>();
   auto *fittedTrack = fittedTrackOwner.get();
   fittedTrack->SetTrackID(track->GetTrackID());
   // Convention here: Kinematics holds the AT-VERTEX (post back-extrapolation)
   // KE — backward-compatible with all existing analysis macros that read
   // `kin.kineticEnergy`. KinematicsXtr holds the AT-FIRST-CLUSTER (pre
   // back-extrap) KE so consumers (e.g. pi_inspect) can show both side by
   // side. θ/φ are the same in both since the smoothed state's angles are
   // at the first cluster and we don't update them during back-extrap.
   fittedTrack->SetKinematics(KE, theta_s, phi_s);
   fittedTrack->SetKinematicsXtr(KE_first_cluster, theta_s, phi_s);
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
