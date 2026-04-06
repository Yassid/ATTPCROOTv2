#include "AtEvent.h"
#include "AtHit.h"
#include "AtHitCluster.h"
#include "AtPatternCircle2D.h"
#include "AtTrack.h"
#include "AtTrackRefiner.h"
#include "AtTrackSeeder.h"
#include "AtTrackTransformer.h"

#include <gtest/gtest.h>

#include <cmath>
#include <initializer_list>
#include <memory>
#include <tuple>
#include <vector>

namespace {

using XYZPoint = ROOT::Math::XYZPoint;

std::shared_ptr<AtHitCluster> MakeCluster(double x, double y, double z, double charge = 1.0, int timeStamp = 0)
{
   auto cluster = std::make_shared<AtHitCluster>();
   cluster->SetPosition({x, y, z});
   cluster->SetCharge(charge);
   cluster->SetTimeStamp(timeStamp);
   return cluster;
}

AtTrack BuildTrack(std::initializer_list<std::tuple<double, double, double>> hits,
                   std::initializer_list<std::tuple<double, double, double>> clusters)
{
   AtTrack track;

   int hitId = 0;
   for (const auto &[x, y, z] : hits) {
      auto hit = std::make_unique<AtHit>();
      hit->SetHitID(hitId++);
      hit->SetPosition({x, y, z});
      hit->SetCharge(1.0);
      hit->SetTimeStamp(hitId);
      track.AddHit(std::move(hit));
   }

   int clusterId = 0;
   for (const auto &[x, y, z] : clusters) {
      auto cluster = MakeCluster(x, y, z, 1.0, clusterId);
      cluster->SetClusterID(clusterId++);
      track.AddClusterHit(cluster);
   }

   return track;
}

AtTrack BuildArcTrack(double radius, XYZPoint center, const std::vector<double> &angles, double zStart, double zStep)
{
   AtTrack track;
   int hitId = 0;
   for (std::size_t i = 0; i < angles.size(); ++i) {
      auto hit = std::make_unique<AtHit>();
      double angle = angles[i];
      hit->SetHitID(hitId++);
      hit->SetPosition({center.X() + radius * std::cos(angle), center.Y() + radius * std::sin(angle), zStart + i * zStep});
      hit->SetCharge(100.0);
      hit->SetTimeStamp(i);
      track.AddHit(std::move(hit));
   }
   return track;
}

void ExpectPointNear(const XYZPoint &lhs, const XYZPoint &rhs, double tol = 1e-9)
{
   EXPECT_NEAR(lhs.X(), rhs.X(), tol);
   EXPECT_NEAR(lhs.Y(), rhs.Y(), tol);
   EXPECT_NEAR(lhs.Z(), rhs.Z(), tol);
}

} // namespace

TEST(AtTrackRefinerTest, OrderClustersAlongTrackStartsAtHighestZAndGreedilyWalksNearestNeighbor)
{
   // Minimal single-track geometry: one outgoing candidate with clusters already
   // lying along a monotonic trajectory in drift Z. This is physically reasonable
   // as a bare ordering fixture, but intentionally omits detector effects.
   AtPATTERN::AtTrackRefiner refiner;
   auto track = BuildTrack(
      {{0.0, 0.0, 100.0}, {0.5, 0.0, 95.0}, {1.0, 0.0, 90.0}, {2.0, 0.0, 80.0}},
      {{1.0, 0.0, 90.0}, {2.0, 0.0, 80.0}, {0.0, 0.0, 100.0}, {0.5, 0.0, 95.0}});

   refiner.OrderClustersAlongTrack(track);

   auto *clusters = track.GetHitClusterArray();
   ASSERT_EQ(clusters->size(), 4u);
   ExpectPointNear(clusters->at(0).GetPosition(), {0.0, 0.0, 100.0});
   ExpectPointNear(clusters->at(1).GetPosition(), {0.5, 0.0, 95.0});
   ExpectPointNear(clusters->at(2).GetPosition(), {1.0, 0.0, 90.0});
   ExpectPointNear(clusters->at(3).GetPosition(), {2.0, 0.0, 80.0});
}

TEST(AtTrackRefinerTest, SelectAndMergeTracksRejectsBeamLikeTracksMergesFragmentsAndDropsIsolatedTracks)
{
   // Synthetic event-topology fixture representing:
   // 1) one primary emerging near the beam axis,
   // 2) one nearby downstream fragment that should merge into that primary,
   // 3) one isolated off-axis fragment that should be rejected,
   // 4) one beam-like straight-through candidate that should fail the angle cut.
   // This is physics-motivated policy coverage, not a detector-response test.
   AtPATTERN::AtTrackRefiner refiner;
   AtTools::AtTrackTransformer transformer;

   AtTrack primary = BuildTrack(
      {{0.0, 0.0, 100.0}, {10.0, 0.0, 85.0}},
      {{0.0, 0.0, 100.0}, {10.0, 0.0, 85.0}});

   AtTrack fragment = BuildTrack(
      {{30.0, 0.0, 82.0}, {12.0, 0.0, 82.0}},
      {{30.0, 0.0, 82.0}, {12.0, 0.0, 82.0}});

   AtTrack isolated = BuildTrack(
      {{120.0, 0.0, 98.0}, {130.0, 0.0, 88.0}},
      {{120.0, 0.0, 98.0}, {130.0, 0.0, 88.0}});

   AtTrack beamLike = BuildTrack(
      {{0.0, 0.0, 110.0}, {0.0, 0.0, 70.0}},
      {{0.0, 0.0, 110.0}, {0.0, 0.0, 70.0}});

   std::vector<AtTrack> tracks;
   tracks.push_back(primary);
   tracks.push_back(fragment);
   tracks.push_back(isolated);
   tracks.push_back(beamLike);

   refiner.SelectAndMergeTracks(tracks, transformer, 5.0, 8.0, 5.0, 6.0, 10.0);

   ASSERT_EQ(tracks.size(), 1u);
   auto &merged = tracks.front();
   EXPECT_EQ(merged.GetHitArray().size(), 4u);
   EXPECT_GE(merged.GetHitClusterArray()->size(), 2u);
   EXPECT_GT(merged.GetHitClusterArray()->front().GetPosition().Z(), merged.GetHitClusterArray()->back().GetPosition().Z());
}

TEST(AtTrackSeederTest, SetTrackInitialParametersCharacterizesCurrentCurvedTrackSeeding)
{
   // Stylized curved-track fixture: a clean circular arc in the transverse plane
   // with monotonic drift Z, representing the idealized trajectory of a charged
   // particle in field. It is suitable for geometric seeding checks, but does not
   // model diffusion, inefficiency, or irregular hit sampling.
   AtPATTERN::AtTrackSeeder seeder;
   auto track = BuildArcTrack(100.0, {0.0, 0.0, 0.0}, {2.30, 2.10, 1.90, 1.70, 1.50, 1.30, 1.10, 0.90}, 20.0, 8.0);

   seeder.SetTrackInitialParameters(track, 0.5, 10, 25);

   auto center = track.GetGeoCenter();
   EXPECT_NEAR(center.first, 0.0, 1.0);
   EXPECT_NEAR(center.second, 0.0, 1.0);
   EXPECT_NEAR(track.GetGeoRadius(), 100.0, 1.0);
   EXPECT_TRUE(std::isnan(track.GetGeoTheta()));
   EXPECT_NE(track.GetGeoPhi(), 0.0);

   auto *pattern = dynamic_cast<const AtPatterns::AtPatternCircle2D *>(track.GetPattern());
   ASSERT_NE(pattern, nullptr);
   EXPECT_NEAR(pattern->GetCenter().X(), 0.0, 1.0);
   EXPECT_NEAR(pattern->GetCenter().Y(), 0.0, 1.0);
   EXPECT_NEAR(pattern->GetRadius(), 100.0, 1.0);
}
