#include "AtHit.h"
#include "AtHitCluster.h"
#include "AtTrack.h"
#include "AtTrackTransformer.h"

#include <gtest/gtest.h>

#include <cmath>
#include <memory>

namespace {

AtTrack BuildTrack(std::initializer_list<std::tuple<double, double, double, double, int>> hits)
{
   AtTrack track;
   for (const auto &[x, y, z, q, t] : hits) {
      auto hit = std::make_unique<AtHit>();
      hit->SetPosition({x, y, z});
      hit->SetCharge(q);
      hit->SetTimeStamp(t);
      track.AddHit(std::move(hit));
   }
   return track;
}

void ExpectNear(double lhs, double rhs, double tol = 1e-9)
{
   EXPECT_NEAR(lhs, rhs, tol);
}

AtHit::XYZVector GetDefaultHitVariance(const AtHit &hit)
{
   constexpr double coefT = 0.00009;
   constexpr double coefL = 0.0000009;
   constexpr double driftVel = 1.0;
   constexpr double tbTime = 0.320;
   constexpr double padResXY = 2.3;
   constexpr double padResZ = padResXY * 1.5;

   double driftTime = hit.GetPosition().Z() / (10.0 * driftVel);
   double varT = 100.0 * coefT * 2.0 * driftTime;
   double varL = 100.0 * coefL * 2.0 * driftTime;
   double tbRes_mm = driftVel * tbTime * 10.0;
   double varTB = tbRes_mm * tbRes_mm / 12.0;
   return {padResXY * padResXY + varT, padResXY * padResXY + varT, padResZ * padResZ + varTB + varL};
}

} // namespace

TEST(AtTrackTransformerCovarianceTest, Smooth3DCurrentBehaviorIsCharacterized)
{
   auto track = BuildTrack({
      {0.0, 0.0, 10.0, 1.0, 0},
      {1.0, 0.0, 20.0, 2.0, 1},
      {20.0, 0.0, 30.0, 3.0, 2},
   });

   AtTools::AtTrackTransformer transformer;
   transformer.ClusterizeSmooth3D(track, 5.0, 10.0);

   auto *clusters = track.GetHitClusterArray();
   ASSERT_EQ(clusters->size(), 2u);

   {
      const auto &cl = clusters->at(0);
      const auto &cov = cl.GetCovMatrix();
      ExpectNear(cl.GetPosition().X(), 0.0);
      ExpectNear(cl.GetPosition().Y(), 0.0);
      ExpectNear(cl.GetPosition().Z(), 10.0);
      ExpectNear(cov(0, 0), 5.308);
      ExpectNear(cov(0, 1), 0.0);
      ExpectNear(cov(0, 2), 0.0);
      ExpectNear(cov(1, 1), 5.308);
      ExpectNear(cov(1, 2), 0.0);
      ExpectNear(cov(2, 2), 12.756013333333334);
   }

   {
      const auto &cl = clusters->at(1);
      const auto &cov = cl.GetCovMatrix();
      ExpectNear(cl.GetPosition().X(), 1.0);
      ExpectNear(cl.GetPosition().Y(), 0.0);
      ExpectNear(cl.GetPosition().Z(), 20.0);
      ExpectNear(cov(0, 0), 5.326);
      ExpectNear(cov(0, 1), 0.0);
      ExpectNear(cov(0, 2), 0.0);
      ExpectNear(cov(1, 1), 5.326);
      ExpectNear(cov(1, 2), 0.0);
      ExpectNear(cov(2, 2), 12.756193333333334);
   }
}

TEST(AtTrackTransformerCovarianceTest, ClusterizeByGroupCurrentBehaviorIsCharacterized)
{
   auto track = BuildTrack({
      {0.0, 0.0, 10.0, 1.0, 0},
      {1.0, 0.0, 20.0, 2.0, 1},
      {2.0, 0.0, 30.0, 3.0, 2},
   });

   AtTools::AtTrackTransformer transformer;
   transformer.ClusterizeByGroup(track, 3);

   auto *clusters = track.GetHitClusterArray();
   ASSERT_EQ(clusters->size(), 1u);

   const auto &cl = clusters->front();
   const auto &cov = cl.GetCovMatrix();
   ExpectNear(cl.GetPosition().X(), 1.3333333333333333);
   ExpectNear(cl.GetPosition().Y(), 0.0);
   ExpectNear(cl.GetPosition().Z(), 20.0);
   ExpectNear(cov(0, 0), 2.0752222222222225);
   ExpectNear(cov(0, 1), 0.0);
   ExpectNear(cov(0, 2), 0.0);
   ExpectNear(cov(1, 1), 2.0752222222222225);
   ExpectNear(cov(1, 2), 0.0);
   ExpectNear(cov(2, 2), 4.9607818518518515);
}

TEST(AtTrackTransformerCovarianceTest, HitClusterOnlineUsesAtHitClusterCovariance)
{
   auto track = BuildTrack({
      {0.0, 0.0, 10.0, 1.0, 0},
      {1.0, 0.0, 20.0, 2.0, 1},
      {2.0, 0.0, 30.0, 3.0, 2},
   });

   AtTools::AtTrackTransformer transformer;
   transformer.SetCovarianceMode(AtTools::AtTrackTransformer::CovarianceMode::HitClusterOnline);
   transformer.ClusterizeByGroup(track, 3);

   auto *clusters = track.GetHitClusterArray();
   ASSERT_EQ(clusters->size(), 1u);

   const auto &cl = clusters->front();
   const auto &cov = cl.GetCovMatrix();

   ExpectNear(cl.GetPosition().X(), 1.3333333333333333);
   ExpectNear(cl.GetPosition().Y(), 0.0);
   ExpectNear(cl.GetPosition().Z(), 20.0);

   AtHitCluster expected;
   for (const auto &hitPtr : track.GetHitArray()) {
      AtHit hit(*hitPtr);
      hit.SetPositionVariance(GetDefaultHitVariance(hit));
      expected.AddHit(hit);
   }

   const auto &expectedCov = expected.GetCovMatrix();
   for (int r = 0; r < 3; ++r) {
      for (int c = 0; c < 3; ++c)
         ExpectNear(cov(r, c), expectedCov(r, c));
   }

   EXPECT_GT(std::abs(cov(0, 0) - 2.0752222222222225), 1e-6);
}
