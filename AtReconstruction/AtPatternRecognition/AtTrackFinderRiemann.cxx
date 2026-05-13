#include "AtTrackFinderRiemann.h"

#include "AtEvent.h"
#include "AtHit.h"
#include "AtPatternCircle2D.h"
#include "AtPatternEvent.h"
#include "AtTrack.h"
#include "AtTrackTransformer.h"

#include <Math/Point3D.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <random>
#include <set>
#include <utility>
#include <vector>

namespace AtPATTERN {

constexpr auto cRED = "\033[1;31m";
constexpr auto cGREEN = "\033[1;32m";
constexpr auto cNORMAL = "\033[0m";

namespace {

struct Circle {
   double cx{0}, cy{0}, r{0};
   bool valid{false};
};

struct ZLine {
   double a{0}, b{0};
   bool valid{false};
};

/// Analytic circle through 3 (x,y) points.
Circle circleThroughThree(const std::array<double, 2> &a, const std::array<double, 2> &b,
                          const std::array<double, 2> &c)
{
   const double ax = a[0], ay = a[1];
   const double bx = b[0], by = b[1];
   const double cx = c[0], cy = c[1];

   const double d = 2.0 * (ax * (by - cy) + bx * (cy - ay) + cx * (ay - by));
   if (std::abs(d) < 1e-9)
      return {};

   const double a2 = ax * ax + ay * ay;
   const double b2 = bx * bx + by * by;
   const double c2 = cx * cx + cy * cy;
   const double ux = (a2 * (by - cy) + b2 * (cy - ay) + c2 * (ay - by)) / d;
   const double uy = (a2 * (cx - bx) + b2 * (ax - cx) + c2 * (bx - ax)) / d;
   const double r = std::hypot(ax - ux, ay - uy);

   if (!(r > 0.5) || !(r < 1.0e4))
      return {};

   return {ux, uy, r, true};
}

inline double pointCircleDist(const Circle &c, double x, double y)
{
   return std::abs(std::hypot(x - c.cx, y - c.cy) - c.r);
}

/// Z = a + b*phi line through the two seeds with the largest |Δphi|.
/// If the phi span is too small (seeds are essentially co-azimuthal), the slope is
/// underdetermined; return invalid so callers know to skip the z residual check.
ZLine zLineFromSeeds(double phi0, double z0, double phi1, double z1, double phi2, double z2)
{
   const double dp01 = std::abs(phi1 - phi0);
   const double dp02 = std::abs(phi2 - phi0);
   const double dp12 = std::abs(phi2 - phi1);
   double pa, za, pb, zb, dp;
   if (dp01 >= dp02 && dp01 >= dp12) {
      pa = phi0; za = z0; pb = phi1; zb = z1; dp = dp01;
   } else if (dp02 >= dp12) {
      pa = phi0; za = z0; pb = phi2; zb = z2; dp = dp02;
   } else {
      pa = phi1; za = z1; pb = phi2; zb = z2; dp = dp12;
   }
   if (dp < 0.05) // seeds span <~3° in phi: slope is too uncertain to filter on
      return {};
   const double b = (zb - za) / (pb - pa);
   const double a = za - b * pa;
   return {a, b, true};
}

/// LSQ refit of z = a + b*phi on a hit set, sorted phi around (cx, cy).
ZLine zLineLSQ(const std::vector<std::array<double, 2>> &xy, const std::vector<double> &z,
                const std::vector<size_t> &idx, const Circle &c)
{
   if (idx.size() < 2)
      return {};
   double sphi = 0, sz = 0, sphi2 = 0, sphiz = 0;
   const double n = static_cast<double>(idx.size());
   for (size_t i : idx) {
      const double phi = std::atan2(xy[i][1] - c.cy, xy[i][0] - c.cx);
      sphi += phi;
      sz += z[i];
      sphi2 += phi * phi;
      sphiz += phi * z[i];
   }
   const double det = n * sphi2 - sphi * sphi;
   if (std::abs(det) < 1e-9)
      return {0.0, 0.0, false}; // phis collapsed; no slope information
   const double b = (n * sphiz - sphi * sz) / det;
   const double a = (sz - b * sphi) / n;
   return {a, b, true};
}

/// Kasa linear-LSQ circle refit on selected indices (Cramer's rule on 3x3 normals).
Circle kasaRefit(const std::vector<std::array<double, 2>> &xy, const std::vector<size_t> &idx)
{
   if (idx.size() < 3)
      return {};
   double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0, sxr = 0, syr = 0, sr = 0;
   for (size_t k : idx) {
      const double x = xy[k][0], y = xy[k][1];
      const double r2 = x * x + y * y;
      sx += x; sy += y;
      sxx += x * x; syy += y * y; sxy += x * y;
      sxr += x * r2; syr += y * r2; sr += r2;
   }
   const double n = static_cast<double>(idx.size());
   const double det = sxx * (syy * n - sy * sy) - sxy * (sxy * n - sy * sx) + sx * (sxy * sy - syy * sx);
   if (std::abs(det) < 1e-9)
      return {};
   const double dA = sxr * (syy * n - sy * sy) - sxy * (syr * n - sy * sr) + sx * (syr * sy - syy * sr);
   const double dB = sxx * (syr * n - sy * sr) - sxr * (sxy * n - sy * sx) + sx * (sxy * sr - syr * sx);
   const double dC = sxx * (syy * sr - syr * sy) - sxy * (sxy * sr - syr * sx) + sxr * (sxy * sy - syy * sx);
   const double A = dA / det, B = dB / det, C = dC / det;
   const double cx = A / 2.0, cy = B / 2.0;
   const double r2 = C + cx * cx + cy * cy;
   if (!(r2 > 0))
      return {};
   return {cx, cy, std::sqrt(r2), true};
}

/// Arc-walk extension: starting from the RANSAC inlier set, walk outward in phi
/// (around the global refined-circle center) from both ends, using a sliding-window
/// LOCAL Kasa refit. Accepts a candidate hit if its perpendicular distance to the
/// local circle is below xyTol AND |z - (a + b*phi_local)| is below zTol. Allows up
/// to maxMiss consecutive failures before stopping a walk direction — bridges small
/// gaps (e.g. dead pads) without latching onto another track.
std::vector<size_t> arcWalkExtend(const std::vector<std::array<double, 2>> &xy,
                                   const std::vector<double> &zVals,
                                   const std::vector<size_t> &inliers,
                                   const std::vector<size_t> &liveAll,
                                   const Circle &globalCircle,
                                   double xyTol, double zTol,
                                   int windowSize = 10, int maxMiss = 5)
{
   if (inliers.size() < (size_t)std::max(3, windowSize / 2))
      return inliers;

   auto phiOf = [&](size_t i) {
      return std::atan2(xy[i][1] - globalCircle.cy, xy[i][0] - globalCircle.cx);
   };

   std::set<size_t> inlierSet(inliers.begin(), inliers.end());
   std::vector<std::pair<double, size_t>> phiInl;
   phiInl.reserve(inliers.size());
   for (size_t i : inliers) phiInl.emplace_back(phiOf(i), i);
   std::sort(phiInl.begin(), phiInl.end());

   std::vector<std::pair<double, size_t>> phiCand;
   phiCand.reserve(liveAll.size());
   for (size_t i : liveAll) {
      if (inlierSet.count(i)) continue;
      phiCand.emplace_back(phiOf(i), i);
   }
   std::sort(phiCand.begin(), phiCand.end());

   auto checkAccept = [&](const std::vector<size_t> &win, size_t cand) -> bool {
      const Circle lc = kasaRefit(xy, win);
      if (!lc.valid) return false;
      if (pointCircleDist(lc, xy[cand][0], xy[cand][1]) >= xyTol) return false;
      const ZLine zl = zLineLSQ(xy, zVals, win, lc);
      if (zl.valid) {
         const double phi = std::atan2(xy[cand][1] - lc.cy, xy[cand][0] - lc.cx);
         if (std::abs(zVals[cand] - (zl.a + zl.b * phi)) >= zTol) return false;
      }
      return true;
   };

   // Forward walk (+phi direction).
   {
      auto it = std::upper_bound(phiCand.begin(), phiCand.end(), phiInl.back(),
                                  [](const std::pair<double, size_t> &a,
                                     const std::pair<double, size_t> &b) {
                                     return a.first < b.first;
                                  });
      int miss = 0;
      for (; it != phiCand.end() && miss <= maxMiss; ++it) {
         const size_t k = std::min((size_t)windowSize, phiInl.size());
         std::vector<size_t> win;
         win.reserve(k);
         for (auto rit = phiInl.end() - k; rit != phiInl.end(); ++rit) win.push_back(rit->second);
         if (checkAccept(win, it->second)) {
            phiInl.emplace_back(it->first, it->second);
            inlierSet.insert(it->second);
            miss = 0;
         } else {
            ++miss;
         }
      }
   }

   // Backward walk (-phi direction).
   {
      auto it = std::lower_bound(phiCand.begin(), phiCand.end(), phiInl.front(),
                                  [](const std::pair<double, size_t> &a,
                                     const std::pair<double, size_t> &b) {
                                     return a.first < b.first;
                                  });
      int miss = 0;
      while (it != phiCand.begin() && miss <= maxMiss) {
         --it;
         const size_t k = std::min((size_t)windowSize, phiInl.size());
         std::vector<size_t> win;
         win.reserve(k);
         for (auto fit = phiInl.begin(); fit != phiInl.begin() + k; ++fit) win.push_back(fit->second);
         if (checkAccept(win, it->second)) {
            phiInl.insert(phiInl.begin(), {it->first, it->second});
            inlierSet.insert(it->second);
            miss = 0;
         } else {
            ++miss;
         }
      }
   }

   std::vector<size_t> out;
   out.reserve(phiInl.size());
   for (auto &p : phiInl) out.push_back(p.second);
   return out;
}

/// Largest contiguous arc around (cx, cy): consecutive sorted phi gaps must stay below maxGap.
/// Returns the surviving subset of input indices (the biggest run on one side of the largest gap).
std::vector<size_t> largestContiguousArc(const std::vector<std::array<double, 2>> &xy,
                                          const std::vector<size_t> &inliers, const Circle &c, double maxGapDeg)
{
   if (inliers.size() < 2)
      return inliers;

   const double maxGap = maxGapDeg * M_PI / 180.0;
   std::vector<std::pair<double, size_t>> phiHits;
   phiHits.reserve(inliers.size());
   for (size_t i : inliers) {
      const double phi = std::atan2(xy[i][1] - c.cy, xy[i][0] - c.cx);
      phiHits.emplace_back(phi, i);
   }
   std::sort(phiHits.begin(), phiHits.end());

   const size_t n = phiHits.size();
   size_t gapStart = 0;
   double maxObservedGap = 0.0;
   for (size_t i = 0; i < n; ++i) {
      const double next = (i + 1 < n) ? phiHits[i + 1].first : (phiHits[0].first + 2 * M_PI);
      const double g = next - phiHits[i].first;
      if (g > maxObservedGap) {
         maxObservedGap = g;
         gapStart = i;
      }
   }
   if (maxObservedGap <= maxGap) {
      std::vector<size_t> out;
      out.reserve(inliers.size());
      for (auto &p : phiHits)
         out.push_back(p.second);
      return out;
   }

   std::vector<size_t> bestRun;
   for (size_t s = 0; s < n; ++s) {
      const size_t start = (gapStart + 1 + s) % n;
      std::vector<size_t> run{phiHits[start].second};
      for (size_t step = 1; step < n; ++step) {
         const size_t cur = (start + step) % n;
         const size_t prev = (start + step - 1) % n;
         double g = phiHits[cur].first - phiHits[prev].first;
         if (cur < prev)
            g += 2 * M_PI;
         if (g > maxGap)
            break;
         run.push_back(phiHits[cur].second);
      }
      if (run.size() > bestRun.size())
         bestRun = std::move(run);
   }
   return bestRun;
}

/// Brute-force 3D kNN: for each hit, return up to k nearest-neighbor indices (excluding self).
std::vector<std::vector<size_t>>
precomputeKNN(const std::vector<std::array<double, 2>> &xy, const std::vector<double> &z, int k)
{
   const size_t N = xy.size();
   std::vector<std::vector<size_t>> nb(N);
   for (size_t i = 0; i < N; ++i) {
      std::vector<std::pair<double, size_t>> dists;
      dists.reserve(N - 1);
      for (size_t j = 0; j < N; ++j) {
         if (j == i)
            continue;
         const double dx = xy[j][0] - xy[i][0];
         const double dy = xy[j][1] - xy[i][1];
         const double dz = z[j] - z[i];
         dists.emplace_back(dx * dx + dy * dy + dz * dz, j);
      }
      const size_t kk = std::min(static_cast<size_t>(k), dists.size());
      if (kk == 0)
         continue;
      std::nth_element(dists.begin(), dists.begin() + kk, dists.end());
      nb[i].reserve(kk);
      for (size_t m = 0; m < kk; ++m)
         nb[i].push_back(dists[m].second);
   }
   return nb;
}

} // namespace

AtTrackFinderRiemann::AtTrackFinderRiemann() : AtPRA() {}

std::unique_ptr<AtPatternEvent> AtTrackFinderRiemann::FindTracks(AtEvent &event)
{
   auto retEvent = std::make_unique<AtPatternEvent>();

   const Int_t nHits = event.GetNumHits();
   if (nHits < fMinHitsPerTrack)
      return retEvent;

   std::vector<std::array<double, 2>> xy;
   std::vector<double> zVals;
   xy.reserve(nHits);
   zVals.reserve(nHits);
   for (Int_t i = 0; i < nHits; ++i) {
      const auto &p = event.GetHit(i).GetPosition();
      xy.push_back({p.X(), p.Y()});
      zVals.push_back(p.Z());
   }

   // Pre-compute k nearest neighbors in 3D for every hit.
   const auto neighbors = precomputeKNN(xy, zVals, fK_NN);

   std::vector<bool> assigned(nHits, false);
   std::vector<AtTrack> tracks;
   std::mt19937 rng(fSeed);

   for (int trk = 0; trk < fMaxTracks; ++trk) {
      std::vector<size_t> live;
      for (size_t i = 0; i < xy.size(); ++i)
         if (!assigned[i])
            live.push_back(i);
      if (static_cast<int>(live.size()) < fMinHitsPerTrack)
         break;

      Circle bestCircle;
      ZLine bestLine;
      std::vector<size_t> bestInliers;

      std::uniform_int_distribution<size_t> pickLive(0, live.size() - 1);

      for (int it = 0; it < fMaxIterations; ++it) {
         const size_t s0 = live[pickLive(rng)];

         // Filter s0's neighbors to those still unassigned.
         std::vector<size_t> liveNeighbors;
         liveNeighbors.reserve(neighbors[s0].size());
         for (size_t nb : neighbors[s0])
            if (!assigned[nb])
               liveNeighbors.push_back(nb);
         if (liveNeighbors.size() < 2)
            continue;

         std::uniform_int_distribution<size_t> pickNb(0, liveNeighbors.size() - 1);
         const size_t s1 = liveNeighbors[pickNb(rng)];
         const size_t s2 = liveNeighbors[pickNb(rng)];
         if (s1 == s2 || s0 == s1 || s0 == s2)
            continue;

         const Circle cand = circleThroughThree(xy[s0], xy[s1], xy[s2]);
         if (!cand.valid)
            continue;

         const double phi0 = std::atan2(xy[s0][1] - cand.cy, xy[s0][0] - cand.cx);
         const double phi1 = std::atan2(xy[s1][1] - cand.cy, xy[s1][0] - cand.cx);
         const double phi2 = std::atan2(xy[s2][1] - cand.cy, xy[s2][0] - cand.cx);
         const ZLine zline = zLineFromSeeds(phi0, zVals[s0], phi1, zVals[s1], phi2, zVals[s2]);

         std::vector<size_t> inl;
         inl.reserve(live.size());
         for (size_t i : live) {
            if (pointCircleDist(cand, xy[i][0], xy[i][1]) >= fInlierDist)
               continue;
            if (zline.valid) {
               const double phi = std::atan2(xy[i][1] - cand.cy, xy[i][0] - cand.cx);
               if (std::abs(zVals[i] - (zline.a + zline.b * phi)) >= fZInlierDist)
                  continue;
            }
            inl.push_back(i);
         }

         if (inl.size() > bestInliers.size()) {
            bestInliers = std::move(inl);
            bestCircle = cand;
            bestLine = zline;
         }
      }

      if (static_cast<int>(bestInliers.size()) < fMinHitsPerTrack)
         break; // No good seed produced a track this pass; stop searching.

      // Snapshot to drop on failure (keeps RANSAC from re-picking the same bad cluster).
      const std::vector<size_t> ransacCluster = bestInliers;
      auto rejectAndContinue = [&]() {
         for (size_t i : ransacCluster)
            assigned[i] = true;
      };

      // Kasa LSQ refit of the xy circle on inliers, plus an LSQ refit of the z line.
      const Circle refinedCircle = kasaRefit(xy, bestInliers);
      if (!refinedCircle.valid) {
         rejectAndContinue();
         continue;
      }
      const ZLine refinedLine = zLineLSQ(xy, zVals, bestInliers, refinedCircle);

      // Re-collect inliers in 3D under the refined circle + refined z line.
      std::vector<size_t> finalInliers;
      finalInliers.reserve(bestInliers.size());
      for (size_t i : live) {
         if (pointCircleDist(refinedCircle, xy[i][0], xy[i][1]) >= fInlierDist)
            continue;
         if (refinedLine.valid) {
            const double phi = std::atan2(xy[i][1] - refinedCircle.cy, xy[i][0] - refinedCircle.cx);
            if (std::abs(zVals[i] - (refinedLine.a + refinedLine.b * phi)) >= fZInlierDist)
               continue;
         }
         finalInliers.push_back(i);
      }
      if (static_cast<int>(finalInliers.size()) < fMinHitsPerTrack) {
         rejectAndContinue();
         continue;
      }
      bestCircle = refinedCircle;
      bestInliers.swap(finalInliers);

      // Single post-check: largest contiguous arc around the fitted center. Smart seeds + 3D
      // inlier check usually produce a single physical arc, but in rare cases the xy inlier
      // expansion picks up hits from another track that shares the (xy) circle but lies in
      // an opposite-quadrant arc; this filter splits those.
      bestInliers = largestContiguousArc(xy, bestInliers, bestCircle, fMaxPhiGap);
      if (static_cast<int>(bestInliers.size()) < fMinHitsPerTrack) {
         rejectAndContinue();
         continue;
      }

      // Arc-walk tail recovery: re-fit locally on a sliding window and walk
      // outward in phi from both ends. Picks up hits where the real trajectory
      // diverges from the global circle (MS, dE/dx, gentle drift) — but stays
      // anchored to the local fit so it doesn't latch onto another track.
      if (fUseArcWalkExtend) {
         std::vector<size_t> liveAll;
         liveAll.reserve(live.size());
         for (size_t i : live) liveAll.push_back(i);
         bestInliers = arcWalkExtend(xy, zVals, bestInliers, liveAll, bestCircle,
                                      fInlierDist, fZInlierDist,
                                      fArcWalkWindow, fArcWalkMaxMiss);
      }

      AtTrack track;
      for (size_t i : bestInliers) {
         track.AddHit(event.GetHit(i));
         assigned[i] = true;
      }
      track.SetTrackID(tracks.size());

      if (fUseArcWalk)
         fTrackTransformer->ClusterizeArcWalk(track, fTargetClusters, fMinHitsPerCluster, fArcWalkKNN);
      else
         fTrackTransformer->ClusterizeSmooth3D(track, fClusterRadius, fClusterDistance);
      OrderClustersAlongTrack(track);

      if (kSetPrunning)
         PruneTrack(track);

      tracks.push_back(std::move(track));
   }

   std::cout << cRED << " AtTrackFinderRiemann (smart-seed circle RANSAC): tracks found " << tracks.size() << cNORMAL
             << "\n";

   if (fUseSelectAndMerge) {
      SelectAndMergeTracks(tracks);
      std::cout << cGREEN << " After selection: " << tracks.size() << " tracks" << cNORMAL << "\n";
   }

   for (auto &t : tracks) {
      const auto &hits = t.GetHitArray();
      if (hits.size() < 3)
         continue;
      std::vector<ROOT::Math::XYZPoint> xyPoints;
      xyPoints.reserve(hits.size());
      for (const auto &h : hits) {
         const auto &p = h->GetPosition();
         xyPoints.emplace_back(p.X(), p.Y(), 0.0);
      }
      auto circle = std::make_unique<AtPatterns::AtPatternCircle2D>();
      circle->DefinePattern(xyPoints);
      t.SetPattern(std::move(circle));
   }

   for (auto &t : tracks)
      if (t.GetHitArray().size() > 0)
         SetTrackInitialParameters(t);

   for (size_t i = 0; i < xy.size(); ++i)
      if (!assigned[i])
         retEvent->AddNoise(event.GetHit(i));

   for (auto &t : tracks)
      retEvent->AddTrack(std::move(t));

   return retEvent;
}

} // namespace AtPATTERN
