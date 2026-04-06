#include "AtPRA.h"

#include "AtHit.h"
#include "AtHitCluster.h"
#include "AtPatternCircle2D.h"
#include "AtPatternEvent.h"
#include "AtTrack.h"

#include <gtest/gtest.h>

#include <initializer_list>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

namespace {

using XYZPoint = ROOT::Math::XYZPoint;

class TestPRA : public AtPATTERN::AtPRA {
public:
   std::unique_ptr<AtPatternEvent> FindTracks(AtEvent &) override { return nullptr; }
   void RunSetTrackInitialParameters(AtTrack &track) { SetTrackInitialParameters(track); }
};

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

std::vector<XYZPoint> GetClusterPositions(AtTrack &track)
{
   std::vector<XYZPoint> positions;
   for (const auto &cluster : *track.GetHitClusterArray())
      positions.push_back(cluster.GetPosition());
   return positions;
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

TEST(AtPRABehaviorTest, OrderClustersAlongTrackStartsAtHighestZAndGreedilyWalksNearestNeighbor)
{
   TestPRA pra;
   auto track = BuildTrack(
      {{0.0, 0.0, 100.0}, {0.5, 0.0, 95.0}, {1.0, 0.0, 90.0}, {2.0, 0.0, 80.0}},
      {{1.0, 0.0, 90.0}, {2.0, 0.0, 80.0}, {0.0, 0.0, 100.0}, {0.5, 0.0, 95.0}});

   pra.OrderClustersAlongTrack(track);

   auto ordered = GetClusterPositions(track);
   ASSERT_EQ(ordered.size(), 4u);
   ExpectPointNear(ordered[0], {0.0, 0.0, 100.0});
   ExpectPointNear(ordered[1], {0.5, 0.0, 95.0});
   ExpectPointNear(ordered[2], {1.0, 0.0, 90.0});
   ExpectPointNear(ordered[3], {2.0, 0.0, 80.0});
}

TEST(AtPRABehaviorTest, SelectAndMergeTracksRejectsBeamLikeTracksMergesFragmentsAndDropsIsolatedTracks)
{
   TestPRA pra;
   pra.SetClusterRadius(5.0);
   pra.SetClusterDistance(8.0);

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

   pra.SelectAndMergeTracks(tracks, 5.0, 6.0, 10.0);

   ASSERT_EQ(tracks.size(), 1u);

   auto &merged = tracks.front();
   EXPECT_EQ(merged.GetHitArray().size(), 4u);
   EXPECT_GE(merged.GetHitClusterArray()->size(), 2u);

   auto *clusters = merged.GetHitClusterArray();
   ASSERT_FALSE(clusters->empty());
   auto first = clusters->front().GetPosition();
   auto last = clusters->back().GetPosition();
   EXPECT_GT(first.Z(), last.Z());
   EXPECT_LT(std::sqrt(first.X() * first.X() + first.Y() * first.Y()), 5.0);
}

TEST(AtPRABehaviorTest, SetTrackInitialParametersCharacterizesCurrentCurvedTrackSeeding)
{
   TestPRA pra;
   auto track = BuildArcTrack(100.0, {0.0, 0.0, 0.0}, {2.30, 2.10, 1.90, 1.70, 1.50, 1.30, 1.10, 0.90}, 20.0, 8.0);

   pra.RunSetTrackInitialParameters(track);

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
