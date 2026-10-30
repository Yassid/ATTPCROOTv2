#include "AtHelix3D.h"

#include "AtHit.h"
#include "AtTrack.h"

#include <Math/Point3D.h>

#include <gtest/gtest.h>

#include <algorithm>
#include <cmath>
#include <random>

using XYZPoint = ROOT::Math::XYZPoint;

namespace {
constexpr double kB = 2.85;    // T
constexpr double kZpad = 1000; // z_lab = 1000 - z_digi

/// A helix launched from the beam axis at z_lab = zv, polar angle thetaDeg, radius R [mm], turning
/// COUNTER-CLOCKWISE in x-y as it travels (the vertex end then shows rotation sense +1, the a1975
/// convention). The track is a band 3 hits wide (+-2 mm radially), one row every 3 mm of arc, swept
/// over `sweep` rad. Charge rises x2 over the last third (a Bragg rise at the stopping end).
/// Hits are stored in the DIGI frame (z_digi = 1000 - z_lab), as the reconstruction produces them.
AtTrack makeHelix(double R, double thetaDeg, double zv, double sweep, unsigned seed = 1, bool shuffle = false)
{
   AtTrack tr;
   std::vector<AtHit> hits;
   const double th = thetaDeg * M_PI / 180, dzdphi = R / std::tan(th);
   const double ds = 3.0 / R; // 3 mm of arc per row
   int id = 0;
   for (double dphi = 0; dphi <= sweep; dphi += ds) {
      const double phi = M_PI + dphi; // centre at (R,0): the start point is the origin
      for (double dr : {-2.0, 0.0, 2.0}) {
         double x = R + (R + dr) * std::cos(phi), y = (R + dr) * std::sin(phi);
         double zLab = zv + dzdphi * dphi;
         double q = 100 * (dphi > 2 * sweep / 3 ? 2.0 : 1.0);
         hits.emplace_back(id++, XYZPoint(x, y, kZpad - zLab), q);
      }
   }
   if (shuffle)
      std::shuffle(hits.begin(), hits.end(), std::mt19937(seed));
   for (auto &h : hits)
      tr.AddHit(h);
   return tr;
}

double expectedBrho(double R, double thetaDeg)
{
   return kB * R / 1000.0 / std::sin(thetaDeg * M_PI / 180);
}
} // namespace

TEST(AtHelix3DTest, RecoversHelixParameters)
{
   AtTools::AtHelix3D h3;
   h3.SetBField(kB);
   auto tr = makeHelix(80, 70, 400, 1.2);
   auto r = h3.Estimate(tr);
   ASSERT_TRUE(r.valid) << "failCode " << r.failCode;
   EXPECT_NEAR(r.brho, expectedBrho(80, 70), 0.03 * expectedBrho(80, 70));
   EXPECT_NEAR(r.theta * 180 / M_PI, 70, 2.0);
   EXPECT_NEAR(r.zLabVertexEnd, 400, 5.0); // the vertex end is the start of the helix
   EXPECT_FALSE(r.backward);
   EXPECT_FALSE(r.ambiguous);
   EXPECT_LT(r.dca, 3.0); // the helix passes through the beam axis
}

TEST(AtHelix3DTest, InputOrderDoesNotMatter)
{
   // the pad grid makes exactly tied 3D distances; ties are broken by POSITION, so any permutation of
   // the hit array gives the same answer up to floating-point summation order
   AtTools::AtHelix3D h3;
   auto a = h3.Estimate(makeHelix(80, 70, 400, 1.2));
   ASSERT_TRUE(a.valid);
   for (unsigned seed : {3u, 7u, 11u, 19u}) {
      auto b = h3.Estimate(makeHelix(80, 70, 400, 1.2, seed, true));
      ASSERT_TRUE(b.valid);
      EXPECT_NEAR(a.brho, b.brho, 1e-6 * a.brho);
      EXPECT_NEAR(a.theta, b.theta, 1e-6);
      EXPECT_NEAR(a.sqrtdEdx, b.sqrtdEdx, 1e-6 * a.sqrtdEdx);
      EXPECT_EQ(a.vertexEnd >= 0, b.vertexEnd >= 0);
      EXPECT_EQ(a.backward, b.backward);
   }
}

TEST(AtHelix3DTest, BackwardTrack)
{
   AtTools::AtHelix3D h3;
   auto r = h3.Estimate(makeHelix(80, 120, 600, 1.2));
   ASSERT_TRUE(r.valid) << "failCode " << r.failCode;
   EXPECT_TRUE(r.backward);
   EXPECT_NEAR(r.theta * 180 / M_PI, 120, 2.0);
   EXPECT_NEAR(r.brho, expectedBrho(80, 120), 0.03 * expectedBrho(80, 120));
}

TEST(AtHelix3DTest, StrayBlobNearVertexIsRemoved)
{
   // a compact spray 20-30 mm inside the circle near the vertex, joined to the track by a gap: it
   // must not become the backbone end nor pull the vertex-end circle small
   AtTools::AtHelix3D h3;
   auto tr = makeHelix(80, 80, 400, 1.2);
   int id = 10000;
   for (int i = 0; i < 12; ++i)
      tr.AddHit(AtHit(id++, XYZPoint(25 + (i % 4) * 1.5, 8 + (i / 4) * 1.5, kZpad - 402), 80));
   auto r = h3.Estimate(tr);
   ASSERT_TRUE(r.valid) << "failCode " << r.failCode;
   EXPECT_GT(r.nDropped, 0);
   EXPECT_NEAR(r.brho, expectedBrho(80, 80), 0.05 * expectedBrho(80, 80));
}

TEST(AtHelix3DTest, TooFewHits)
{
   AtTools::AtHelix3D h3;
   h3.SetMinHits(30);
   auto r = h3.Estimate(makeHelix(80, 70, 400, 0.05));
   EXPECT_FALSE(r.valid);
   EXPECT_EQ(r.failCode, 1);
}
