#include "AtSpyralPID.h"

#include "AtHit.h"
#include "AtHitCluster.h"
#include "AtTrack.h"

#include <algorithm>
#include <cmath>
#include <numeric>

using namespace AtTools;

// ---------------------------------------------------------------------------
//  Smoothing spline  (scipy make_smoothing_spline equivalent, Green & Silverman)
// ---------------------------------------------------------------------------

bool AtSmoothingSpline::Fit(const std::vector<double> &x, const std::vector<double> &y, double lam)
{
   fValid = false;
   const int n = (int)x.size();
   if (n < 5 || (int)y.size() != n)
      return false;
   // strictly increasing knots required (else scipy raises -> Spyral skips event)
   for (int i = 0; i + 1 < n; ++i)
      if (x[i + 1] <= x[i])
         return false;

   const int m = n - 2; // interior knots
   std::vector<double> h(n - 1);
   for (int i = 0; i + 1 < n; ++i)
      h[i] = x[i + 1] - x[i];

   // Q (n x m), R (m x m, tridiagonal)
   std::vector<double> Q((size_t)n * m, 0.0);
   std::vector<double> R((size_t)m * m, 0.0);
   for (int k = 1; k <= n - 2; ++k) { // interior knot index
      int c = k - 1;
      Q[(size_t)(k - 1) * m + c] += 1.0 / h[k - 1];
      Q[(size_t)k * m + c] += -(1.0 / h[k - 1] + 1.0 / h[k]);
      Q[(size_t)(k + 1) * m + c] += 1.0 / h[k];
      R[(size_t)c * m + c] = (h[k - 1] + h[k]) / 3.0;
   }
   for (int a = 0; a + 1 < m; ++a) {
      double v = h[a + 1] / 6.0;
      R[(size_t)a * m + (a + 1)] = v;
      R[(size_t)(a + 1) * m + a] = v;
   }

   // M = R + lam Q^T Q is symmetric pentadiagonal (Q columns are 3-banded), so build
   // and solve it within bandwidth 2 -- O(m) instead of O(m^3). Q column i is nonzero
   // only at rows {i, i+1, i+2}.
   std::vector<double> M((size_t)m * m, 0.0);
   std::vector<double> Qty(m, 0.0);
   auto Qel = [&](int r, int c) { return Q[(size_t)r * m + c]; };
   for (int i = 0; i < m; ++i) {
      for (int j = i; j <= i + 2 && j < m; ++j) {
         double s = 0.0;
         for (int r = j; r <= i + 2; ++r) // overlap of nonzero rows of cols i and j
            s += Qel(r, i) * Qel(r, j);
         double val = R[(size_t)i * m + j] + lam * s;
         M[(size_t)i * m + j] = val;
         M[(size_t)j * m + i] = val;
      }
      double sy = 0.0;
      for (int r = i; r <= i + 2; ++r)
         sy += Qel(r, i) * y[r];
      Qty[i] = sy;
   }

   // banded LU (bandwidth 2, no pivoting -- M is SPD pentadiagonal)
   std::vector<double> gamma(m, 0.0);
   for (int k = 0; k < m; ++k) {
      double piv = M[(size_t)k * m + k];
      if (std::fabs(piv) < 1e-15)
         return false;
      for (int i = k + 1; i <= k + 2 && i < m; ++i) {
         double fct = M[(size_t)i * m + k] / piv;
         if (fct == 0.0)
            continue;
         for (int j = k; j <= k + 2 && j < m; ++j)
            M[(size_t)i * m + j] -= fct * M[(size_t)k * m + j];
         Qty[i] -= fct * Qty[k];
      }
   }
   for (int i = m - 1; i >= 0; --i) {
      double s = Qty[i];
      for (int j = i + 1; j <= i + 2 && j < m; ++j)
         s -= M[(size_t)i * m + j] * gamma[j];
      gamma[i] = s / M[(size_t)i * m + i];
   }

   // g = y - lam Q gamma  (fitted/smoothed values at the knots)
   fG.assign(n, 0.0);
   for (int r = 0; r < n; ++r) {
      double s = 0.0;
      for (int i = 0; i < m; ++i)
         s += Q[(size_t)r * m + i] * gamma[i];
      fG[r] = y[r] - lam * s;
   }
   // second derivatives at knots (natural spline: ends zero, interior = gamma)
   fY2.assign(n, 0.0);
   for (int k = 1; k <= n - 2; ++k)
      fY2[k] = gamma[k - 1];
   fX = x;
   fValid = true;
   return true;
}

double AtSmoothingSpline::Eval(double t) const
{
   const int n = (int)fX.size();
   if (n == 0)
      return 0.0;
   if (t <= fX.front())
      t = fX.front();
   if (t >= fX.back())
      t = fX.back();
   int lo = 0, hi = n - 1;
   while (hi - lo > 1) {
      int mid = (lo + hi) / 2;
      if (fX[mid] > t)
         hi = mid;
      else
         lo = mid;
   }
   double h = fX[hi] - fX[lo];
   if (h == 0.0)
      return fG[lo];
   double A = (fX[hi] - t) / h;
   double B = (t - fX[lo]) / h;
   return A * fG[lo] + B * fG[hi] + ((A * A * A - A) * fY2[lo] + (B * B * B - B) * fY2[hi]) * (h * h) / 6.0;
}

// ---------------------------------------------------------------------------
//  Geometry helpers
// ---------------------------------------------------------------------------

bool AtSpyralPID::LeastSquaresCircle(const std::vector<double> &x, const std::vector<double> &y, double &xc, double &yc,
                                     double &radius)
{
   const int n = (int)x.size();
   if (n < 3)
      return false;
   double mx = std::accumulate(x.begin(), x.end(), 0.0) / n;
   double my = std::accumulate(y.begin(), y.end(), 0.0) / n;
   double Suv = 0, Suu = 0, Svv = 0, Suuv = 0, Suvv = 0, Suuu = 0, Svvv = 0;
   for (int i = 0; i < n; ++i) {
      double u = x[i] - mx, v = y[i] - my;
      Suv += u * v;
      Suu += u * u;
      Svv += v * v;
      Suuv += u * u * v;
      Suvv += u * v * v;
      Suuu += u * u * u;
      Svvv += v * v * v;
   }
   double det = Suu * Svv - Suv * Suv;
   if (std::fabs(det) < 1e-15)
      return false;
   double bu = (Suuu + Suvv) * 0.5, bv = (Suuv + Svvv) * 0.5;
   double uc = (bu * Svv - bv * Suv) / det;
   double vc = (Suu * bv - Suv * bu) / det;
   xc = uc + mx;
   yc = vc + my;
   double meanR = 0.0;
   for (int i = 0; i < n; ++i)
      meanR += std::sqrt((x[i] - xc) * (x[i] - xc) + (y[i] - yc) * (y[i] - yc));
   radius = meanR / n;
   return true;
}

int AtSpyralPID::Direction(const std::vector<double> &xin, const std::vector<double> &yin,
                           const std::vector<double> &z) const
{
   const int n = (int)xin.size();
   if (n < 2)
      return -1;
   std::vector<double> x = xin, y = yin;
   if (n > 5) {
      AtSmoothingSpline sx, sy;
      if (sx.Fit(z, xin, 10.0) && sy.Fit(z, yin, 10.0)) {
         x = sx.FittedValues();
         y = sy.FittedValues();
      }
   }
   double cx, cy, r;
   if (!LeastSquaresCircle(x, y, cx, cy, r))
      return -1;
   std::vector<double> ang(n);
   for (int i = 0; i < n; ++i) {
      double a = std::atan2(y[i] - cy, x[i] - cx);
      if (a < 0)
         a += 2 * M_PI;
      ang[i] = a * 180.0 / M_PI;
   }
   int pos = 0, tot = n - 1;
   for (int i = 0; i + 1 < n; ++i)
      if (ang[i + 1] - ang[i] > 0)
         ++pos;
   double posFrac = tot > 0 ? double(pos) / tot : 0.0;
   double negFrac = 1.0 - posFrac;
   if (posFrac > fDirThreshold)
      return 0; // FORWARD
   if (negFrac > fDirThreshold)
      return 1; // BACKWARD
   return -1;   // NONE
}

// ---------------------------------------------------------------------------
//  Estimator  (faithful port of estimate_physics / estimate_physics_pass)
// ---------------------------------------------------------------------------

AtSpyralResult AtSpyralPID::Estimate(AtTrack &track) const
{
   AtSpyralResult res;

   // Stamp the identity and quality fields FIRST, before any of the early returns below, so that
   // even a track whose PID fails still carries them. A gate is drawn on the whole landscape,
   // including the rejects, and having these here means that can be done from the PID output
   // alone -- no need to reopen the pattern file to recover which track this was or how big it was.
   res.trackID = track.GetTrackID();
   {
      auto *hc = track.GetHitClusterArray();
      res.nClusters = hc ? static_cast<int>(hc->size()) : 0;
   }

   // gather trajectory points (x,y,z,charge), sorted ascending in z.
   // Faithful to Spyral: use the AtHit point cloud (GetHitArray). Optionally fall
   // back to the coarser merged hit-clusters.
   struct P {
      double x, y, z, q;
   };
   std::vector<P> pts;
   if (fUseHits) {
      auto &hits = track.GetHitArray();
      if (hits.empty())
         return res;
      pts.reserve(hits.size());
      for (auto &h : hits) {
         if (!h)
            continue;
         auto p = h->GetPosition();
         // Spyral integrates the pulse; AtHit::GetCharge() is the peak amplitude. Using the
         // trace integral matches Spyral's dE/dx scale (~x4). Opt-in via SetUseTraceIntegral.
         double q = (fUseTraceIntegral && h->GetTraceIntegral() > 0) ? h->GetTraceIntegral() : h->GetCharge();
         pts.push_back({p.X(), p.Y(), p.Z(), q});
      }
   } else {
      auto *clusters = track.GetHitClusterArray();
      if (!clusters || clusters->empty())
         return res;
      pts.reserve(clusters->size());
      for (auto &c : *clusters) {
         auto p = c.GetPosition();
         pts.push_back({p.X(), p.Y(), p.Z(), c.GetCharge()});
      }
   }
   std::sort(pts.begin(), pts.end(), [](const P &a, const P &b) { return a.z < b.z; });
   // Collapse exact z-ties: AT-TPC z is discretized into time buckets, so several
   // clusters can share a z. The smoothing spline needs strictly increasing z, so
   // merge tied points into one (charge-weighted x,y; summed charge). This is the
   // point-cloud preprocessing that lets Spyral's spline-based estimator run here.
   {
      std::vector<P> merged;
      for (size_t i = 0; i < pts.size();) {
         size_t j = i;
         double wx = 0, wy = 0, qs = 0;
         while (j < pts.size() && pts[j].z == pts[i].z) {
            double w = pts[j].q > 0 ? pts[j].q : 1.0;
            wx += w * pts[j].x;
            wy += w * pts[j].y;
            qs += pts[j].q;
            ++j;
         }
         double wsum = 0;
         for (size_t k = i; k < j; ++k)
            wsum += (pts[k].q > 0 ? pts[k].q : 1.0);
         merged.push_back({wx / wsum, wy / wsum, pts[i].z, qs});
         i = j;
      }
      pts.swap(merged);
   }
   int n = (int)pts.size();
   if (n < fMinPoints) // min_total_trajectory_points, on the usable (collapsed) cloud
      { res.failCode = 1; return res; }

   std::vector<double> x(n), y(n), z(n), q(n);
   for (int i = 0; i < n; ++i) {
      x[i] = pts[i].x;
      y[i] = pts[i].y;
      z[i] = pts[i].z;
      q[i] = pts[i].q;
   }

   // direction from clustering-style winding analysis
   int direction = Direction(x, y, z);
   if (direction < 0)
      { res.failCode = 2; return res; }

   // smoothing splines (lam=fSmoothing) on x(z), y(z), charge(z); replace data
   AtSmoothingSpline xSpline, ySpline, cSpline;
   if (!xSpline.Fit(z, x, fSmoothing) || !ySpline.Fit(z, y, fSmoothing) || !cSpline.Fit(z, q, fSmoothing))
      { res.failCode = 3; return res; } // spline failure (e.g. multivalued z) -> skip, like Spyral
   x = xSpline.FittedValues();
   y = ySpline.FittedValues();
   q = cSpline.FittedValues();

   // flip if backward so the trajectory starts at the vertex end
   if (direction == 1) {
      std::reverse(x.begin(), x.end());
      std::reverse(y.begin(), y.end());
      std::reverse(z.begin(), z.end());
      std::reverse(q.begin(), q.end());
   }

   // vertex guess = first point
   double vx = x[0], vy = y[0], vz = 0.0;

   // first arc: furthest point from vertex in rho, over points 1..n-1
   auto argmaxRho = [&](int start) {
      int best = start;
      double bestv = -1;
      for (int i = start; i < n; ++i) {
         double rr = std::hypot(x[i] - vx, y[i] - vy);
         if (rr > bestv) {
            bestv = rr;
            best = i;
         }
      }
      return best; // returns the actual point index of the maximum
   };
   // Spyral: rho over cluster_data[1:], argmax j, first_arc = data[:j+1] (0..j),
   // i.e. excludes the furthest point by one. Replicate: maxIdx is the position
   // in the [1:] slice -> first_arc end index = that position (0-based) = j.
   {
      int bestj = 0;
      double bestv = -1;
      for (int i = 1; i < n; ++i) {
         double rr = std::hypot(x[i] - vx, y[i] - vy);
         if (rr > bestv) {
            bestv = rr;
            bestj = i - 1; // position within data[1:]
         }
      }
      int maximum = bestj; // first_arc = indices 0..maximum
      res.dbgMax1 = bestj;
      // fit circle to first arc
      std::vector<double> fx(x.begin(), x.begin() + maximum + 1), fy(y.begin(), y.begin() + maximum + 1);
      double cx, cy, radius;
      if (!LeastSquaresCircle(fx, fy, cx, cy, radius))
         { res.failCode = 4; return res; }

      // re-estimate vertex: circle point nearest the beam axis (origin in xy)
      double dc = std::hypot(cx, cy);
      if (dc > 1e-9) {
         vx = cx * (1.0 - radius / dc);
         vy = cy * (1.0 - radius / dc);
      }

      // recompute rho to refined vertex over ALL points; new first arc
      std::vector<double> rho(n);
      int maxAll = 0;
      double maxv = -1;
      for (int i = 0; i < n; ++i) {
         rho[i] = std::hypot(x[i] - vx, y[i] - vy);
         if (rho[i] > maxv) {
            maxv = rho[i];
            maxAll = i;
         }
      }
      int maximum2 = maxAll; // first_arc = 0..maximum2
      res.dbgMax2 = maxAll;

      // linear fit rho vs z over the first 50% of the first arc (>=10 pts)
      int testIndex = std::max(10, (int)(maximum2 * 0.5));
      if (testIndex > n)
         testIndex = n;
      res.dbgTestIdx = testIndex;
      if (testIndex < 2)
         { res.failCode = 5; return res; }
      double sz = 0, sr = 0;
      for (int i = 0; i < testIndex; ++i) {
         sz += z[i];
         sr += rho[i];
      }
      double mz = sz / testIndex, mr = sr / testIndex;
      double Sxy = 0, Sxx = 0;
      for (int i = 0; i < testIndex; ++i) {
         Sxy += (z[i] - mz) * (rho[i] - mr);
         Sxx += (z[i] - mz) * (z[i] - mz);
      }
      if (Sxx == 0.0)
         { res.failCode = 6; return res; }
      double slope = Sxy / Sxx;
      res.dbgSlope = slope;
      double intercept = mr - slope * mz;
      if (slope == 0.0)
         { res.failCode = 7; return res; }
      double vertexRho = std::hypot(vx, vy);
      vz = -intercept / slope;

      if (vertexRho > fBeamRegionRadius)
         { res.failCode = 8; return res; }

      double polar = std::atan(slope);
      if ((polar > 0.0 && direction == 1) || (polar < 0.0 && direction == 0))
         { res.failCode = 9; return res; } // direction/polar inconsistent
      if (direction == 1)
         polar += M_PI;

      double azimuthal = std::atan2(vy - cy, vx - cx);
      if (azimuthal < 0)
         azimuthal += 2.0 * M_PI;
      azimuthal += M_PI * 0.5;
      if (azimuthal > 2.0 * M_PI)
         azimuthal -= 2.0 * M_PI;

      double brho = fBField * radius * 0.001 / std::fabs(std::sin(polar));
      if (std::isnan(brho))
         brho = 0.0;

      // integrate charge over the first arc until leaving the small-pad region
      double chargeDep = q[0];
      int cutoff = -1;
      for (int idx = 0; idx < maximum2; ++idx) {
         if (std::hypot(x[idx + 1], y[idx + 1]) > fSmallPadRadius) {
            cutoff = idx + 1;
            break;
         }
         chargeDep += q[idx + 1];
      }
      if (chargeDep == q[0])
         { res.failCode = 10; return res; }

      int cutIdx = (cutoff == -1) ? maximum2 : cutoff; // Spyral: -1 -> last of first_arc
      // arc length: 1000-point spline line integral over the first-arc z-range
      const int NP = 1000;
      double z0 = z[0], z1 = z[cutIdx];
      double arclength = 0.0;
      double px = xSpline.Eval(z0), py = ySpline.Eval(z0), pz = z0;
      for (int k = 1; k < NP; ++k) {
         double zk = z0 + (z1 - z0) * k / (NP - 1);
         double xk = xSpline.Eval(zk), yk = ySpline.Eval(zk);
         arclength += std::sqrt((xk - px) * (xk - px) + (yk - py) * (yk - py) + (zk - pz) * (zk - pz));
         px = xk;
         py = yk;
         pz = zk;
      }
      if (arclength <= 0.0)
         { res.failCode = 11; return res; }

      double dEdx = chargeDep / arclength;

      res.vertex = ROOT::Math::XYZPoint(vx, vy, vz);
      res.center = ROOT::Math::XYZPoint(cx, cy, vz);
      res.polar = polar;
      res.azimuthal = azimuthal;
      res.radius = radius;
      res.brho = brho;
      res.dE = chargeDep;
      res.arclength = arclength;
      res.dEdx = dEdx;
      res.sqrtdEdx = std::sqrt(std::fabs(dEdx));
      res.direction = direction;
      res.nPoints = n;
      res.valid = true;
   }
   return res;
}
