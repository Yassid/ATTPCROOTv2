#include "AtELossCATIMA.h"
#include "AtHitCluster.h"
#include "AtKinematics.h"
#include "AtPropagator.h"
#include "AtTrack.h"
#include "AtTrackTransformer.h"

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <TMath.h>
#include <TMatrixD.h>

#include <TRandom.h>

#include <cmath>
#include <gtest/gtest.h>
#include <iomanip>
#include <iostream>
#include <memory>
#include <tuple>
#include <vector>

#include "OpenKF/kalman_filter/TrackFitterUKF.h"

using ROOT::Math::Polar3DVector;
using ROOT::Math::XYZPoint;
using ROOT::Math::XYZVector;

// ===========================================================================
// Generate a realistic proton track using the propagator, then smear the
// hit positions to simulate digitization effects.
// ===========================================================================
class ClusteringScanTest : public testing::Test {
protected:
   static constexpr double kMass = 938.272;
   static constexpr double kCharge = 1.602176634e-19;
   static constexpr double kBz = 2.85;
   static constexpr double kGasDensity = 3.553e-5;

   // True initial state
   XYZPoint kVertex{0, 0, 170};                         // mm
   XYZVector kTrueMom{9.35463, -45.4279, 8.26042};      // MeV/c
   double kTrueP = kTrueMom.R();                         // ~47.1 MeV/c

   // Generate raw hits along a propagated track with position smearing.
   // hitSpacing_mm controls how far apart the hits are placed along the track
   // to simulate the density of pad hits from digitization (~1-2 mm).
   std::vector<XYZPoint> GenerateHits(double hitSpacing_mm = 1.5, double smear_mm = 1.0)
   {
      auto eloss = std::make_unique<AtTools::AtELossCATIMA>(kGasDensity);
      eloss->SetProjectile(1, 1, 1);
      std::vector<std::tuple<int, int, int>> mat;
      mat.push_back({1, 1, 1});
      eloss->SetMaterial(mat);

      AtTools::AtPropagator prop(kCharge, kMass, std::move(eloss));
      prop.SetBField({0, 0, kBz});
      prop.SetState(kVertex, kTrueMom);

      AtTools::AtRK4Stepper stepper;
      std::vector<XYZPoint> hits;
      hits.push_back(kVertex);

      double accumulated = 0;
      for (int i = 0; i < 10000; i++) { // Max iterations
         prop.PropagateOneStep(stepper);
         auto state = prop.GetState();
         if (AtTools::Kinematics::KE(state.fMom, kMass) < 0.05)
            break;

         accumulated += (state.fPos - state.fLastPos).R();
         if (accumulated < hitSpacing_mm)
            continue;
         accumulated = 0;

         // Smear position to simulate pad + diffusion resolution
         double sx = gRandom->Gaus(0, smear_mm);
         double sy = gRandom->Gaus(0, smear_mm);
         double sz = gRandom->Gaus(0, smear_mm * 0.5);
         XYZPoint smeared(state.fPos.X() + sx, state.fPos.Y() + sy, state.fPos.Z() + sz);
         hits.push_back(smeared);
      }
      return hits;
   }

   // Build an AtTrack from raw hit positions
   AtTrack BuildTrack(const std::vector<XYZPoint> &hits)
   {
      AtTrack track;
      for (size_t i = 0; i < hits.size(); i++) {
         auto hit = std::make_unique<AtHit>();
         hit->SetPosition(hits[i]);
         hit->SetCharge(100.0); // Arbitrary charge
         hit->SetTimeStamp(i);
         track.AddHit(std::move(hit));
      }
      // Set geometry params for Brho seed
      track.SetGeoRadius(350.0);
      track.SetGeoTheta(kTrueMom.Theta());
      track.SetGeoPhi(kTrueMom.Phi());
      return track;
   }

   // Create UKF instance
   std::unique_ptr<kf::TrackFitterUKF> CreateUKF()
   {
      auto eloss = std::make_unique<AtTools::AtELossCATIMA>(kGasDensity);
      eloss->SetProjectile(1, 1, 1);
      std::vector<std::tuple<int, int, int>> mat;
      mat.push_back({1, 1, 1});
      eloss->SetMaterial(mat);

      AtTools::AtPropagator propagator(kCharge, kMass, std::move(eloss));
      propagator.SetBField({0, 0, kBz});
      auto stepper = std::make_unique<AtTools::AtRK4AdaptiveStepper>();
      auto ukf = std::make_unique<kf::TrackFitterUKF>(std::move(propagator), std::move(stepper));
      ukf->setParameters(1e-3, 2.0, 0.0);
      ukf->fEnableEnStraggling = true;
      return ukf;
   }

   // Run UKF on hit clusters and return (momentum_error%, converged)
   struct FitResult {
      double momErr;
      double meanResid;
      int nClusters;
      bool converged;
   };

   FitResult RunUKF(kf::TrackFitterUKF &ukf, const std::vector<AtHitCluster> &clusters)
   {
      FitResult result{0, 0, static_cast<int>(clusters.size()), false};
      if (clusters.size() < 5)
         return result;

      // Get ordered cluster positions (they're already sorted by ClusterizeSmooth3D)
      std::vector<XYZPoint> pts;
      for (auto &cl : clusters)
         pts.push_back(cl.GetPosition());

      // Filter by minimum spacing
      std::vector<XYZPoint> filtered;
      filtered.push_back(pts[0]);
      for (size_t i = 1; i < pts.size(); i++) {
         if ((pts[i] - filtered.back()).R() >= 2.0) // 2mm min spacing
            filtered.push_back(pts[i]);
      }
      result.nClusters = filtered.size();
      if (result.nClusters < 5)
         return result;

      // Setup
      double sigma_pos = 2.0;
      double sigma_mom = 0.15 * kTrueP;
      double sigma_ang = 5.0 * M_PI / 180.0;
      TMatrixD cov(6, 6);
      cov.Zero();
      for (int i = 0; i < 3; i++)
         cov(i, i) = sigma_pos * sigma_pos;
      cov(3, 3) = sigma_mom * sigma_mom;
      cov(4, 4) = sigma_ang * sigma_ang;
      cov(5, 5) = sigma_ang * sigma_ang;

      TMatrixD measCov(3, 3);
      measCov.Zero();
      for (int i = 0; i < 3; i++)
         measCov(i, i) = sigma_pos * sigma_pos;

      ukf.SetInitialState(filtered[0], kTrueMom, cov);
      ukf.SetMeasCov(measCov);

      try {
         for (size_t i = 1; i < filtered.size(); i++) {
            ukf.predictUKF(filtered[i]);
            ukf.correctUKF(filtered[i]);
         }
         ukf.smoothUKF();
      } catch (...) {
         return result;
      }

      result.converged = true;
      double pReco = ukf.GetSmoothedStates()[0][3];
      result.momErr = (pReco - kTrueP) / kTrueP * 100;

      auto &smoothed = ukf.GetSmoothedStates();
      double sumResid = 0;
      int n = std::min(smoothed.size(), filtered.size()) - 1;
      for (int i = 1; i <= n; i++) {
         XYZPoint sp(smoothed[i][0], smoothed[i][1], smoothed[i][2]);
         sumResid += (sp - filtered[i]).R();
      }
      result.meanResid = n > 0 ? sumResid / n : 0;

      return result;
   }
};

// ===========================================================================
// Test: scan clustering radius and distance parameters
// ===========================================================================
TEST_F(ClusteringScanTest, RadiusDistanceScan)
{
   // Generate raw hits
   auto hits = GenerateHits(1.5, 1.0);
   ASSERT_GT(hits.size(), 50u) << "Not enough hits generated";

   auto ukf = CreateUKF();

   struct Config {
      double radius;
      double distance;
   };
   std::vector<Config> configs = {
      {5.0, 10.0},  {5.0, 15.0},  {10.0, 10.0}, {10.0, 15.0}, {10.0, 20.0},
      {15.0, 15.0}, {15.0, 20.0}, {15.0, 30.0}, {15.0, 30.5}, {20.0, 20.0},
      {20.0, 30.0}, {20.0, 40.0}, {25.0, 30.0}, {25.0, 40.0},
   };

   std::cout << "\n"
             << std::setw(10) << "radius" << std::setw(10) << "distance" << std::setw(10) << "clusters"
             << std::setw(12) << "mom_err(%)" << std::setw(10) << "resid(mm)" << std::setw(10) << "status"
             << std::endl;
   std::cout << std::string(62, '-') << std::endl;

   int nConverged = 0;
   for (auto &cfg : configs) {
      AtTrack track = BuildTrack(hits);
      AtTools::AtTrackTransformer transformer;
      transformer.ClusterizeSmooth3D(track, cfg.radius, cfg.distance);

      auto *clusters = track.GetHitClusterArray();
      auto result = RunUKF(*ukf, *clusters);

      std::string status = result.converged ? "OK" : "FAIL";
      if (result.converged && std::abs(result.momErr) > 5.0)
         status = "WARN";

      std::cout << std::setw(10) << cfg.radius << std::setw(10) << cfg.distance << std::setw(10) << result.nClusters
                << std::setw(12) << (result.converged ? std::to_string(result.momErr).substr(0, 6) : "---")
                << std::setw(10) << (result.converged ? std::to_string(result.meanResid).substr(0, 5) : "---")
                << std::setw(10) << status << std::endl;

      if (result.converged)
         nConverged++;
   }

   // At least half should converge
   EXPECT_GT(nConverged, static_cast<int>(configs.size()) / 2)
      << "Too many clustering configurations failed to converge";
}

// ===========================================================================
// Test: statistical scan — run multiple trials with smeared hits
// ===========================================================================
TEST_F(ClusteringScanTest, StatisticalScan)
{
   auto ukf = CreateUKF();

   struct Config {
      double radius;
      double distance;
      std::string label;
   };
   std::vector<Config> configs = {
      {5.0, 10.0, "r5_d10"},   {10.0, 15.0, "r10_d15"}, {15.0, 20.0, "r15_d20"},
      {15.0, 30.5, "r15_d30"}, {20.0, 30.0, "r20_d30"}, {25.0, 40.0, "r25_d40"},
   };

   int nTrials = 20;

   std::cout << "\n--- Statistical scan (" << nTrials << " trials per config) ---\n" << std::endl;
   std::cout << std::setw(12) << "config" << std::setw(8) << "fit" << std::setw(8) << "fail" << std::setw(12)
             << "bias(%)" << std::setw(12) << "rms(%)" << std::setw(10) << "resid" << std::setw(10) << "clusters"
             << std::endl;
   std::cout << std::string(72, '-') << std::endl;

   for (auto &cfg : configs) {
      std::vector<double> errs;
      std::vector<double> resids;
      std::vector<int> nCls;
      int nFail = 0;

      for (int t = 0; t < nTrials; t++) {
         auto hits = GenerateHits(1.5, 1.0);
         AtTrack track = BuildTrack(hits);
         AtTools::AtTrackTransformer transformer;
         transformer.ClusterizeSmooth3D(track, cfg.radius, cfg.distance);

         auto result = RunUKF(*ukf, *track.GetHitClusterArray());
         if (result.converged && std::abs(result.momErr) < 50) {
            errs.push_back(result.momErr);
            resids.push_back(result.meanResid);
            nCls.push_back(result.nClusters);
         } else {
            nFail++;
         }
      }

      if (errs.empty()) {
         std::cout << std::setw(12) << cfg.label << std::setw(8) << 0 << std::setw(8) << nFail
                   << "  (all failed)" << std::endl;
         continue;
      }

      double mean = 0, rms = 0, mResid = 0, mCl = 0;
      for (auto v : errs)
         mean += v;
      for (auto v : resids)
         mResid += v;
      for (auto v : nCls)
         mCl += v;
      mean /= errs.size();
      mResid /= resids.size();
      mCl /= nCls.size();
      for (auto v : errs)
         rms += (v - mean) * (v - mean);
      rms = std::sqrt(rms / errs.size());

      std::cout << std::setw(12) << cfg.label << std::setw(8) << errs.size() << std::setw(8) << nFail
                << std::setw(12) << std::fixed << std::setprecision(2) << mean << std::setw(12) << rms
                << std::setw(10) << mResid << std::setw(10) << std::setprecision(0) << mCl << std::endl;

      // Basic sanity: bias should be reasonable and convergence rate > 50%
      EXPECT_LT(std::abs(mean), 10.0) << cfg.label << " has excessive bias";
   }
}
