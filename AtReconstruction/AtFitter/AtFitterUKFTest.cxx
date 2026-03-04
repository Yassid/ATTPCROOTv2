#include "AtFitterUKF.h"

#include "AtELossCATIMA.h"
#include "AtFitMetadata.h"
#include "AtFittedTrack.h"
#include "AtHitCluster.h"
#include "AtTrack.h"

#include <Math/Point3D.h>
#include <Math/Vector3D.h>
#include <TMath.h>

#include <cmath>
#include <gtest/gtest.h>
#include <memory>
#include <tuple>
#include <vector>

// Expose protected GetFittedTrack for testing.
class AtFitterUKFTestable : public EventFit::AtFitterUKF {
public:
   using EventFit::AtFitterUKF::AtFitterUKF;
   using EventFit::AtFitterUKF::GetFittedTrack; // make public
};

// ---------------------------------------------------------------------------
// Test fixture
// ---------------------------------------------------------------------------

class AtFitterUKFFixture : public testing::Test {
protected:
   // Proton properties
   static constexpr double kMass_MeV = 938.272;         // MeV/c^2
   static constexpr double kCharge_C = 1.602176634e-19; // Coulombs (Z=1)
   static constexpr double kBz = 2.85;                  // Tesla
   static constexpr double kTrueKE_MeV = 47.0;          // Expected initial KE from simulation

   // Cluster positions (mm) — hardcoded from hits.txt (pointsToCluster=5, rows 0,5,10,...,35)
   // Original values in cm are multiplied by 10 to convert to mm.
   const std::vector<ROOT::Math::XYZPoint> kClusters{
      {-3.40046e-03, -1.49863e-03, 1.0018}, // row 0
      {0.777186, -4.85911, 1.87907},        // row 5
      {1.10508, -9.76677, 2.76864},         // row 10
      {0.982273, -14.6827, 3.66392},        // row 15
      {0.559072, -18.5943, 4.38033},        // row 20
      {-0.40074, -23.4325, 5.18767},        // row 25
      {-1.80375, -28.1617, 5.99325},        // row 30
      {-3.62341, -32.7446, 6.81099},        // row 35
   };

   std::unique_ptr<AtFitterUKFTestable> fFitter;

   // Build a CATIMA energy-loss model for proton in 300 Torr H2.
   // CATIMA is analytic — no external file needed, full energy range supported.
   static std::unique_ptr<AtTools::AtELossCATIMA> MakeELoss()
   {
      auto eloss = std::make_unique<AtTools::AtELossCATIMA>(3.553e-5); // density g/cm^3
      eloss->SetProjectile(1, 1, 1);                                   // proton: A, Z, charge
      std::vector<std::tuple<int, int, int>> mat;
      mat.push_back({1, 1, 1}); // H2: A=1, Z=1, stoich=1
      eloss->SetMaterial(mat);
      return eloss;
   }

   void SetUp() override
   {
      fFitter = std::make_unique<AtFitterUKFTestable>(kCharge_C, kMass_MeV, MakeELoss());
      fFitter->SetBField({0, 0, kBz});
      fFitter->SetUKFParameters(1e-3, 2.0, 0.0);
      fFitter->SetMeasurementSigma(1.0);
      fFitter->SetMinClusters(3);
      // Energy straggling is on by default (AtFitterUKF default); CATIMA supports it.
   }

   // Build an AtTrack populated with the hardcoded clusters.
   // GeoRadius/Theta/Phi are tuned so the Brho seed gives ~300 MeV/c (KE~47 MeV for proton).
   AtTrack BuildTrack(int numClusters = -1) const
   {
      AtTrack track;
      // Circle-fit parameters: radius=350 mm, theta≈π/2, phi=−π/2 (pointing −Y, +Z)
      track.SetGeoRadius(350.0); // mm → brho = 2.85*0.350/sin(1.48) ≈ 1.00 T·m → p≈300 MeV/c
      track.SetGeoTheta(1.48);   // ~84.8°
      track.SetGeoPhi(-TMath::Pi() / 2.0);

      int n = (numClusters < 0) ? static_cast<int>(kClusters.size()) : numClusters;
      for (int i = 0; i < n && i < static_cast<int>(kClusters.size()); ++i) {
         auto cl = std::make_shared<AtHitCluster>();
         cl->SetPosition(kClusters[i]);
         track.AddClusterHit(cl);
      }
      return track;
   }
};

// ---------------------------------------------------------------------------
// Test: construction does not throw
// ---------------------------------------------------------------------------

TEST_F(AtFitterUKFFixture, Construction)
{
   EXPECT_NO_THROW(AtFitterUKFTestable fitter(kCharge_C, kMass_MeV, MakeELoss()));
}

// ---------------------------------------------------------------------------
// Test: tracks with too few clusters are skipped (GetFittedTrack returns nullptr)
// ---------------------------------------------------------------------------

TEST_F(AtFitterUKFFixture, SkipsTrackWithFewClusters)
{
   fFitter->SetMinClusters(5);
   AtTrack track = BuildTrack(2); // only 2 clusters, below threshold
   AtFittedTrack *result = nullptr;
   EXPECT_NO_THROW(result = fFitter->GetFittedTrack(&track));
   EXPECT_EQ(result, nullptr);
}

// ---------------------------------------------------------------------------
// Test: fitted track has physically reasonable kinematics
// ---------------------------------------------------------------------------

TEST_F(AtFitterUKFFixture, FittedTrackHasReasonableKinematics)
{
   AtTrack track = BuildTrack();
   AtFittedTrack *result = nullptr;
   ASSERT_NO_THROW(result = fFitter->GetFittedTrack(&track));
   ASSERT_NE(result, nullptr);

   auto kin = result->GetKinematics();
   EXPECT_GT(kin.kineticEnergy, 0.0);
   EXPECT_GT(kin.theta, 0.0);
   EXPECT_LT(kin.theta, TMath::Pi());

   EXPECT_FALSE(result->GetSmoothedPositions().empty());

   delete result;
}

// ---------------------------------------------------------------------------
// Test: fitted KE is within 20% of the true initial KE (~47 MeV)
// ---------------------------------------------------------------------------

TEST_F(AtFitterUKFFixture, FittedTrackKinematicsWithinTolerance)
{
   AtTrack track = BuildTrack();
   AtFittedTrack *result = nullptr;
   ASSERT_NO_THROW(result = fFitter->GetFittedTrack(&track));
   ASSERT_NE(result, nullptr);

   auto kin = result->GetKinematics();
   double relErr = std::abs(kin.kineticEnergy - kTrueKE_MeV) / kTrueKE_MeV;
   EXPECT_LT(relErr, 0.20) << "Fitted KE = " << kin.kineticEnergy << " MeV, expected ~" << kTrueKE_MeV << " MeV";

   delete result;
}

// ---------------------------------------------------------------------------
// Test: smoothed positions vector has the expected size (nClusters - 1)
// ---------------------------------------------------------------------------

TEST_F(AtFitterUKFFixture, SmoothedPositionsPopulated)
{
   AtTrack track = BuildTrack();
   int nClusters = static_cast<int>(kClusters.size());
   AtFittedTrack *result = nullptr;
   ASSERT_NO_THROW(result = fFitter->GetFittedTrack(&track));
   ASSERT_NE(result, nullptr);

   EXPECT_EQ(static_cast<int>(result->GetSmoothedPositions().size()), nClusters - 1);

   delete result;
}

// ---------------------------------------------------------------------------
// Test: fit metadata reports convergence for a good track
// ---------------------------------------------------------------------------

TEST_F(AtFitterUKFFixture, MetadataConvergedSet)
{
   AtTrack track = BuildTrack();
   AtFittedTrack *result = nullptr;
   ASSERT_NO_THROW(result = fFitter->GetFittedTrack(&track));
   ASSERT_NE(result, nullptr);

   auto &meta = result->GetTrackMetadata();
   ASSERT_NE(meta.get(), nullptr);
   EXPECT_TRUE(meta->GetFitConverged());

   delete result;
}

// ---------------------------------------------------------------------------
// Round-trip tests on AtFittedTrack — no physics needed
// ---------------------------------------------------------------------------

TEST(AtFittedTrackRoundTrip, KinematicsRoundTrip)
{
   AtFittedTrack track;
   track.SetKinematics(42.5, 1.2, -0.7);

   auto kin = track.GetKinematics();
   EXPECT_DOUBLE_EQ(kin.kineticEnergy, 42.5);
   EXPECT_DOUBLE_EQ(kin.theta, 1.2);
   EXPECT_DOUBLE_EQ(kin.phi, -0.7);
}

TEST(AtFittedTrackRoundTrip, SmoothedPositionsRoundTrip)
{
   AtFittedTrack track;
   std::vector<ROOT::Math::XYZPoint> pts{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
   track.SetSmoothedPositions(pts);

   const auto &got = track.GetSmoothedPositions();
   ASSERT_EQ(got.size(), pts.size());
   for (size_t i = 0; i < pts.size(); ++i) {
      EXPECT_DOUBLE_EQ(got[i].X(), pts[i].X());
      EXPECT_DOUBLE_EQ(got[i].Y(), pts[i].Y());
      EXPECT_DOUBLE_EQ(got[i].Z(), pts[i].Z());
   }
}
