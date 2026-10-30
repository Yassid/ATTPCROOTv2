/// @file AtHelix3DCore.h
/// @brief The algorithm behind AtTools::AtHelix3D, header-only so the analysis bench macros
///        (macro/Unpack_HDF5/a1975/UKF/pid/helix3d.h includes this file) and the framework run
///        the SAME code. Developed and validated on a1975 16C+d data, 2026-10-30; see AtHelix3D.h.
///
/// Pipeline: (1) BLOB REMOVAL -- cut the 3D minimum-spanning tree at a multiple of the track's own
/// median edge, drop only small pieces that are NOT on the main piece's circle, capped;
/// (2) BACKBONE -- the tree diameter orders the hits along the track (no z sort, input-order free);
/// (3) CENTRELINE -- charge-weighted centroids in bins of backbone length; (4) END HELICES -- a
/// robust Taubin circle + z(phi) line at EACH end, whole track shrunk to half a turn, with a radius
/// sanity guard and a vertex-on-axis constraint; (5) VERTEX END -- 2-of-3 majority of rotation
/// sense, Spyral direction and a Bragg rise; (6) PID from the vertex end (firstArcPID).
#ifndef ATHELIX3DCORE_H
#define ATHELIX3DCORE_H
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

namespace AtTools {
namespace Helix3D {


struct P {
   double x, y, z, q;
   int id{-1}; // index in the CALLER's hit array (set by estimate, survives blob removal)
};

struct End {
   bool ok{false};
   int nUsed{0};
   double xc{0}, yc{0}, R{0};
   double dphi{0};   // phase swept over the window, signed (travelling inward from this end)
   int sense{0};     // sign of dphi
   double dzdphi{0}; // z_lab per radian
   double dca{1e9};  // closest approach of the circle to the beam axis [mm]
   double gap{1e9};  // arc length from the window start back to that closest approach [mm]
   double vx{0}, vy{0}, vz{0};
   double theta{-1}; // polar of travel inward from this end [rad], z_lab frame
   double brho{0};   // T m
   double qdens{0};  // charge per mm of 3D path over the window
   double qEdge{0};  // charge per mm over the outermost 10 % of hits (Bragg probe)
   double path{0};   // 3D path length of the window [mm]
   double sWin{0};   // backbone length covered by the window, from this end [mm]
   double resCut{1e9}; // robust-fit residual cut [mm]; points farther off the circle are outliers
   int first{0};       // position (in this end's order) of the first INLIER -- the effective end
   bool guarded{false}; // end circle replaced by the whole-track circle (see the sanity guard)
   bool axisConstrained{false}; // end circle refitted through the beam axis (vertex-on-axis constraint)
};

struct Result {
   bool ok{false};
   int n{0}, nBackbone{0};
   double length{0}; // backbone 3D length [mm]
   End e[2];         // e[0] starts at backbone end A, e[1] at end B
   std::vector<int> order; // hit indices, A -> B
   std::vector<double> s;  // backbone coordinate of order[i], from end A [mm]
   int nDropped{0};        // hits removed as small blobs
   std::vector<P> hits;    // the hit set `order` indexes into (lab frame, blobs removed)
   std::vector<P> cen;     // centreline: charge-weighted centroids along the backbone, A -> B (q = summed)
   std::vector<double> sc; // backbone coordinate of each centroid, from end A [mm]
};

/// Taubin algebraic circle fit (Newton solve, after N. Chernov). Much less biased than Kasa toward
/// small radii on short noisy arcs, which is exactly the a1975 near-90-deg case.
inline bool taubin(const std::vector<P> &h, const std::vector<int> &idx, double &xc, double &yc, double &R)
{
   const int n = idx.size();
   if (n < 3) return false;
   double mx = 0, my = 0;
   for (int i : idx) { mx += h[i].x; my += h[i].y; }
   mx /= n; my /= n;
   double Mxx = 0, Myy = 0, Mxy = 0, Mxz = 0, Myz = 0, Mzz = 0;
   for (int i : idx) {
      double x = h[i].x - mx, y = h[i].y - my, z = x * x + y * y;
      Mxy += x * y; Mxx += x * x; Myy += y * y; Mxz += x * z; Myz += y * z; Mzz += z * z;
   }
   Mxx /= n; Myy /= n; Mxy /= n; Mxz /= n; Myz /= n; Mzz /= n;
   double Mz = Mxx + Myy, Cov_xy = Mxx * Myy - Mxy * Mxy, Var_z = Mzz - Mz * Mz;
   double A3 = 4 * Mz, A2 = -3 * Mz * Mz - Mzz;
   double A1 = Var_z * Mz + 4 * Cov_xy * Mz - Mxz * Mxz - Myz * Myz;
   double A0 = Mxz * (Mxz * Myy - Myz * Mxy) + Myz * (Myz * Mxx - Mxz * Mxy) - Var_z * Cov_xy;
   double A22 = A2 + A2, A33 = A3 + A3 + A3, x = 0, y = A0;
   for (int it = 0; it < 99; ++it) {
      double Dy = A1 + x * (A22 + A33 * x), xnew = x - y / Dy;
      if (!std::isfinite(xnew) || xnew == x) break;
      double ynew = A0 + xnew * (A1 + xnew * (A2 + xnew * A3));
      if (std::abs(ynew) >= std::abs(y)) break;
      x = xnew; y = ynew;
   }
   double DET = x * x - x * Mz + Cov_xy;
   if (std::abs(DET) < 1e-12) return false;
   double Xc = (Mxz * (Myy - x) - Myz * Mxy) / DET / 2, Yc = (Myz * (Mxx - x) - Mxz * Mxy) / DET / 2;
   xc = Xc + mx; yc = Yc + my;
   R = std::sqrt(Xc * Xc + Yc * Yc + Mz);
   return std::isfinite(R) && R > 0;
}

/// Robust Taubin: refit without points more than max(cutMin, 3*rms) off the circle, up to 3 times.
/// Stray hit clusters attached to a real track (noise blobs, a neighbour's fragment) otherwise pull a
/// short window's circle to a fraction of the true radius (a1975 (d,p): R 17-36 mm for ~70 mm tracks,
/// 3D Brho 0.44x Spyral on the deuterons the 3D PID misread as protons, 2026-10-30).
inline bool robustCircle(const std::vector<P> &h, std::vector<int> idx, double &xc, double &yc, double &R,
                         std::vector<int> &inl, double &cut, double cutMin = 4.0)
{
   inl = idx;
   for (int it = 0; it < 3; ++it) {
      if (!taubin(h, inl, xc, yc, R)) return false;
      double s2 = 0;
      for (int i : inl) { double d = std::hypot(h[i].x - xc, h[i].y - yc) - R; s2 += d * d; }
      cut = std::max(cutMin, 3 * std::sqrt(s2 / inl.size()));
      std::vector<int> keep;
      for (int i : idx)
         if (std::abs(std::hypot(h[i].x - xc, h[i].y - yc) - R) < cut) keep.push_back(i);
      if (keep.size() < 3) return false;
      if (keep.size() == inl.size()) break;
      inl.swap(keep);
   }
   return taubin(h, inl, xc, yc, R);
}

// MST backbone order. Returns false for degenerate input.
/// Strict weak order on hit POSITION, used to break ties wherever an index would otherwise decide
/// (Prim's start node and equal-length edges, equal backbone positions). It makes the result
/// independent of the order of the input hit array -- the pad grid produces exactly tied distances.
inline bool posLess(const P &a, const P &b)
{
   if (a.x != b.x)
      return a.x < b.x;
   if (a.y != b.y)
      return a.y < b.y;
   if (a.z != b.z)
      return a.z < b.z;
   return a.q < b.q;
}

/// Prim's minimum spanning tree with position-based tie-breaking. Returns parent and squared edge.
inline void primMST(const std::vector<P> &h, std::vector<int> &par, std::vector<double> &best)
{
   const int n = h.size();
   best.assign(n, std::numeric_limits<double>::max());
   par.assign(n, -1);
   std::vector<char> in(n, 0);
   int start = 0;
   for (int i = 1; i < n; ++i)
      if (posLess(h[i], h[start]))
         start = i;
   best[start] = 0;
   for (int it = 0; it < n; ++it) {
      int u = -1;
      for (int i = 0; i < n; ++i)
         if (!in[i] && (u < 0 || best[i] < best[u] || (best[i] == best[u] && posLess(h[i], h[u]))))
            u = i;
      in[u] = 1;
      for (int v = 0; v < n; ++v) {
         if (in[v])
            continue;
         double dx = h[u].x - h[v].x, dy = h[u].y - h[v].y, dz = h[u].z - h[v].z, d = dx * dx + dy * dy + dz * dz;
         if (d < best[v] || (d == best[v] && par[v] >= 0 && posLess(h[u], h[par[v]]))) {
            best[v] = d;
            par[v] = u;
         }
      }
   }
}

/// Label the pieces of the Euclidean MST after cutting edges longer than maxEdge; returns the count.
inline int pieces(const std::vector<P> &h, double maxEdge, std::vector<int> &comp)
{
   const int n = h.size();
   std::vector<double> best;
   std::vector<int> par;
   primMST(h, par, best);
   // maxEdge < 0 means RELATIVE: cut at |maxEdge| x the median MST edge of this track (>= 4 mm).
   // Hit spacing is anisotropic -- a steep track jumps 8-10 mm in z between neighbouring pads, a
   // near-90-deg one has hits 3 mm apart -- so no fixed 3D length separates blobs from track.
   double cutLen = maxEdge;
   if (maxEdge < 0) {
      std::vector<double> e;
      for (int v = 0; v < n; ++v) if (par[v] >= 0) e.push_back(std::sqrt(best[v]));
      std::sort(e.begin(), e.end());
      cutLen = std::max(4.0, -maxEdge * (e.empty() ? 0 : e[e.size() / 2]));
   }
   std::vector<std::vector<int>> adj(n);
   for (int v = 0; v < n; ++v)
      if (par[v] >= 0 && std::sqrt(best[v]) <= cutLen) { adj[v].push_back(par[v]); adj[par[v]].push_back(v); }
   comp.assign(n, -1);
   int nc = 0;
   for (int v = 0; v < n; ++v) {
      if (comp[v] >= 0) continue;
      std::vector<int> st{v};
      comp[v] = nc;
      while (!st.empty()) { int u = st.back(); st.pop_back(); for (int x : adj[u]) if (comp[x] < 0) { comp[x] = nc; st.push_back(x); } }
      ++nc;
   }
   return nc;
}

/// maxEdge: MST edges longer than this are CUT and only the largest connected piece is kept, so a
/// stray cluster joined to the track by a gap cannot become a backbone END (2026-10-30: such blobs
/// set the vertex-end window and halved Brho). Hits outside the kept piece are left out of `order`.
inline bool backbone(const std::vector<P> &h, std::vector<int> &order, int &nBack, double &len,
                     std::vector<double> *sOrd = nullptr, double maxEdge = 1e9)
{
   const int n = h.size();
   if (n < 5) return false;
   std::vector<double> best;
   std::vector<int> par;
   primMST(h, par, best);
   std::vector<std::vector<std::pair<int, double>>> adj(n);
   for (int v = 0; v < n; ++v)
      if (par[v] >= 0) {
         double d = std::sqrt(best[v]);
         if (d > maxEdge) continue; // cut: gap
         adj[v].push_back({par[v], d});
         adj[par[v]].push_back({v, d});
      }
   // largest connected piece
   std::vector<int> comp(n, -1);
   int nComp = 0, bestComp = 0, bestSize = 0;
   for (int v = 0; v < n; ++v) {
      if (comp[v] >= 0) continue;
      int size = 0;
      std::vector<int> st{v};
      comp[v] = nComp;
      while (!st.empty()) {
         int u = st.back(); st.pop_back(); ++size;
         for (auto [x, w] : adj[u]) if (comp[x] < 0) { comp[x] = nComp; st.push_back(x); }
      }
      if (size > bestSize) { bestSize = size; bestComp = nComp; }
      ++nComp;
   }
   if (bestSize < 5) return false;
   int seed = -1; // position-minimal hit of the kept piece (input-order free)
   for (int i = 0; i < n; ++i)
      if (comp[i] == bestComp && (seed < 0 || posLess(h[i], h[seed])))
         seed = i;
   auto far = [&](int s, std::vector<double> &dist, std::vector<int> &from) {
      dist.assign(n, -1); from.assign(n, -1);
      std::vector<int> st{s};
      dist[s] = 0;
      int arg = s;
      while (!st.empty()) {
         int u = st.back(); st.pop_back();
         if (dist[u] > dist[arg] || (dist[u] == dist[arg] && posLess(h[u], h[arg])))
            arg = u;
         for (auto [v, w] : adj[u])
            if (dist[v] < 0) { dist[v] = dist[u] + w; from[v] = u; st.push_back(v); }
      }
      return arg;
   };
   std::vector<double> d1, d2;
   std::vector<int> f1, f2;
   int a = far(seed, d1, f1);
   int b = far(a, d2, f2);
   len = d2[b];
   std::vector<int> path;
   for (int v = b; v >= 0; v = f2[v]) path.push_back(v);
   std::reverse(path.begin(), path.end()); // a -> b
   nBack = path.size();
   // every hit takes the backbone index of the node it hangs from (multi-source flood)
   std::vector<int> root(n, -1);
   std::vector<double> rd(n, 0);
   std::vector<int> st;
   for (int k = 0; k < (int)path.size(); ++k) { root[path[k]] = k; st.push_back(path[k]); }
   while (!st.empty()) {
      int u = st.back(); st.pop_back();
      for (auto [v, w] : adj[u])
         if (root[v] < 0) { root[v] = root[u]; rd[v] = rd[u] + w; st.push_back(v); }
   }
   order.clear();
   for (int i = 0; i < n; ++i) if (comp[i] == bestComp) order.push_back(i);
   std::sort(order.begin(), order.end(), [&](int p, int q) {
      if (root[p] != root[q])
         return root[p] < root[q];
      if (rd[p] != rd[q])
         return rd[p] < rd[q];
      return posLess(h[p], h[q]);
   });
   if (sOrd) {
      sOrd->resize(order.size());
      for (size_t i = 0; i < order.size(); ++i) (*sOrd)[i] = d2[path[root[order[i]]]]; // backbone length from end A
   }
   return true;
}

/// s = backbone arc coordinate of each entry of ord, measured from THIS end (increasing).
/// The window always spans at least minLen mm of track (or the whole track): near the beam axis the
/// hits are dense ACROSS the track, and a window of a fixed number of hits there covers a few mm and
/// fits a meaningless few-mm circle (measured on a1975 (d,p) data, 2026-10-30).
/// rExcl: centroids closer than this to the beam axis do not enter the CIRCLE fit (they still set the
/// window and the phases). The start of an a1975 track, 10-40 mm off axis, is often a diffuse spray
/// (beam-induced charge) that pulls the vertex-end circle to R 14-30 mm (2026-10-30, 19 of the 29
/// remaining 3D proton-band deuterons). Ignored if fewer than 6 window points would remain.
inline End endHelix(const std::vector<P> &h, const std::vector<int> &ord, const std::vector<double> &s, double B,
                    int minPts = 10, double minLen = 40.0, double rExcl = 0.0, double guardFrac = 0.6,
                    double dcaMax = 10.0, double rStartAxis = 45.0)
{
   End E;
   const int n = ord.size();
   if (n < minPts) return E;
   const double total = s.back() - s.front();
   const double needLen = std::min(total, std::max(minLen, 0.25 * total));
   int wMin = n;
   for (int i = 0; i < n; ++i)
      if (s[i] - s[0] >= needLen) { wMin = i + 1; break; }
   wMin = std::max(wMin, minPts);
   // Start from the WHOLE track and shrink to half a turn only if it makes more than that. Growing
   // from a short window fails on wide tracks: a few cm of a 10-20 mm wide band (3-4 pads across,
   // a1975 data near 90 deg) fits a small circle round the blob, the phase sweeps pi quickly round
   // that wrong circle, and the window never grows (2026-10-30: 55/159 Brho collapsed 3-20x).
   int w = n;
   double xc = 0, yc = 0, R = 0;
   std::vector<double> ph;
   std::vector<int> inl;
   double cut = 1e9;
   auto offAxis = [&](std::vector<int> idx) {
      std::vector<int> o;
      for (int i : idx) if (std::hypot(h[i].x, h[i].y) >= rExcl) o.push_back(i);
      return o.size() >= 6 ? o : idx;
   };
   double R0 = 0, sw0 = 0;
   for (int pass = 0; pass < 4; ++pass) {
      std::vector<int> idx = offAxis(std::vector<int>(ord.begin(), ord.begin() + w));
      if (!robustCircle(h, idx, xc, yc, R, inl, cut)) return E;
      if (pass == 0)
         R0 = R; // the WHOLE-track circle (w = n on pass 0)
      // unwrap the phase along the backbone over the WHOLE track, then cut at half a turn
      ph.assign(n, 0);
      double prev = std::atan2(h[ord[0]].y - yc, h[ord[0]].x - xc);
      ph[0] = prev;
      for (int i = 1; i < n; ++i) {
         double a = std::atan2(h[ord[i]].y - yc, h[ord[i]].x - xc), d = a - prev;
         while (d > M_PI) d -= 2 * M_PI;
         while (d < -M_PI) d += 2 * M_PI;
         ph[i] = ph[i - 1] + d;
         prev = a;
      }
      if (pass == 0) sw0 = std::abs(ph[n - 1] - ph[0]);
      int wNew = n;
      for (int i = 0; i < n; ++i)
         if (std::abs(ph[i] - ph[0]) > M_PI) { wNew = i; break; }
      wNew = std::max(wNew, wMin);
      if (wNew == w) break;
      w = wNew;
   }
   std::vector<int> idx = offAxis(std::vector<int>(ord.begin(), ord.begin() + w));
   if (!robustCircle(h, idx, xc, yc, R, inl, cut)) return E;
   // SANITY GUARD: a stopping particle's radius only SHRINKS along its path, so an end circle much
   // smaller than the whole-track circle on a track of under 3/4 turn means the end window was
   // captured by a spray of hits (a1975 vertex region), not by curvature. Use the whole track then.
   if (R0 > 0 && R < guardFrac * R0 && sw0 < 1.5 * M_PI) {
      w = n;
      std::vector<int> all = offAxis(std::vector<int>(ord.begin(), ord.end()));
      if (!robustCircle(h, all, xc, yc, R, inl, cut)) return E;
      E.guarded = true;
   }
   // VERTEX-ON-AXIS CONSTRAINT: if this end STARTS near the beam axis but the free circle misses the
   // axis by more than dcaMax, refit the circle THROUGH the axis (x^2+y^2 = 2ax + 2by, linear in a,b).
   // The vertex is on the axis, so the true vertex-end circle must pass near it. On the 29 remaining
   // 3D proton-band deuterons the free fit missed by 14-30 mm and R_3D was 0.63 of d_max/2, while
   // Spyral's radius was consistent with d_max on 24/29 (fit-free arbiter, 2026-10-30).
   {
      double r0 = 1e9;
      for (int i : inl) r0 = std::min(r0, std::hypot(h[i].x, h[i].y));
      double dcaFree = std::abs(std::hypot(xc, yc) - R);
      if (dcaFree > dcaMax && r0 < rStartAxis && inl.size() >= 3) {
         double Sxx = 0, Sxy = 0, Syy = 0, Sxr = 0, Syr = 0;
         for (int i : inl) {
            double x = h[i].x, y = h[i].y, r2 = x * x + y * y;
            Sxx += x * x; Sxy += x * y; Syy += y * y; Sxr += x * r2; Syr += y * r2;
         }
         double det = Sxx * Syy - Sxy * Sxy;
         if (std::abs(det) > 1e-9) {
            double a = 0.5 * (Sxr * Syy - Syr * Sxy) / det, bb = 0.5 * (Syr * Sxx - Sxr * Sxy) / det;
            xc = a; yc = bb; R = std::hypot(a, bb);
            E.axisConstrained = true;
         }
      }
   }
   // re-unwrap round the final circle, then use INLIERS only (first..last inlier define the window)
   {
      double prev = std::atan2(h[ord[0]].y - yc, h[ord[0]].x - xc);
      ph[0] = prev;
      for (int i = 1; i < n; ++i) {
         double a = std::atan2(h[ord[i]].y - yc, h[ord[i]].x - xc), d = a - prev;
         while (d > M_PI) d -= 2 * M_PI;
         while (d < -M_PI) d += 2 * M_PI;
         ph[i] = ph[i - 1] + d;
         prev = a;
      }
   }
   std::vector<char> isIn(n, 0);
   {
      std::vector<char> mark(h.size(), 0);
      for (int i : inl) mark[i] = 1;
      for (int i = 0; i < w; ++i) isIn[i] = mark[ord[i]];
   }
   int iFirst = -1, iLast = -1, nIn = 0;
   double sp = 0, sz = 0, spp = 0, spz = 0, sq = 0;
   for (int i = 0; i < w; ++i) {
      if (!isIn[i]) continue;
      if (iFirst < 0) iFirst = i;
      iLast = i;
      ++nIn;
      double p = ph[i], z = h[ord[i]].z;
      sp += p; sz += z; spp += p * p; spz += p * z; sq += h[ord[i]].q;
   }
   if (nIn < 3) return E;
   double den = nIn * spp - sp * sp;
   if (std::abs(den) < 1e-12) return E;
   double b = (nIn * spz - sp * sz) / den, a = (sz - b * sp) / nIn;
   E.resCut = cut;
   E.first = iFirst;
   E.dphi = ph[iLast] - ph[iFirst];
   E.sense = E.dphi > 0 ? 1 : -1;
   E.xc = xc; E.yc = yc; E.R = R; E.dzdphi = b; E.nUsed = w; E.sWin = s[iLast] - s[0];
   double dc = std::hypot(xc, yc);
   E.dca = std::abs(dc - R);
   // closest approach on the circle, reached by going BACKWARD (against the inward travel)
   double phv = std::atan2(-yc, -xc);
   double ph0 = ph[iFirst];
   while ((ph0 - phv) * E.sense < 0) phv -= 2 * M_PI * E.sense;
   while ((ph0 - phv) * E.sense >= 2 * M_PI) phv += 2 * M_PI * E.sense;
   E.gap = R * std::abs(ph0 - phv);
   E.vx = xc + R * std::cos(phv); E.vy = yc + R * std::sin(phv); E.vz = a + b * phv;
   // inward travel: dphi/dt has sign `sense`, so dz/dt has sign of b*sense
   double bz = b * E.sense;
   E.theta = std::acos(bz / std::sqrt(bz * bz + R * R));
   double sth = std::sin(E.theta);
   E.brho = sth > 1e-3 ? B * R / 1000.0 / sth : 0;
   E.path = std::sqrt(b * b + R * R) * std::abs(E.dphi);
   E.qdens = E.path > 0 ? sq / E.path : 0;
   // Bragg probe: outermost 10 % of hits by backbone order, charge per mm of chord path
   int k = std::max(3, n / 10);
   double qe = 0, pe = 0;
   for (int i = 0; i < k && i < n; ++i) {
      qe += h[ord[i]].q;
      if (i) {
         auto &u = h[ord[i - 1]], &v = h[ord[i]];
         pe += std::sqrt((u.x - v.x) * (u.x - v.x) + (u.y - v.y) * (u.y - v.y) + (u.z - v.z) * (u.z - v.z));
      }
   }
   E.qEdge = pe > 0 ? qe / pe : 0;
   E.ok = true;
   return E;
}

/// End-fit options (defaults = the adopted configuration of 2026-10-30).
struct EndOpt {
   double rExcl{0.0}, guardFrac{0.6}, dcaMax{10.0}, rStartAxis{45.0}, minLen{40.0};
   double blobLinMin{0.0}; // keep a small piece if its linearity (lambda1/sum) exceeds this (0 = off)
   double blobOffCircle{6.0}; // keep a small piece if it lies within this (mm, mean |d-R|) of the largest
                              // piece's circle -- chopped-off TRACK, not a blob (0 = off)
};

/// Linearity of a point set: largest eigenvalue of its 3D covariance over the trace (1 = a line).
/// Closed-form eigenvalues of the symmetric 3x3 matrix (trigonometric method), so AtTools needs no
/// ROOT::Matrix dependency.
inline double linearity(const std::vector<P> &h, const std::vector<int> &idx)
{
   const int n = idx.size();
   if (n < 3)
      return 0;
   double m[3] = {0, 0, 0};
   for (int i : idx) {
      m[0] += h[i].x;
      m[1] += h[i].y;
      m[2] += h[i].z;
   }
   for (double &v : m)
      v /= n;
   double c[3][3] = {{0}};
   for (int i : idx) {
      double d[3] = {h[i].x - m[0], h[i].y - m[1], h[i].z - m[2]};
      for (int a = 0; a < 3; ++a)
         for (int b = 0; b < 3; ++b)
            c[a][b] += d[a] * d[b] / n;
   }
   const double tr = c[0][0] + c[1][1] + c[2][2];
   if (tr <= 0)
      return 0;
   const double p1 = c[0][1] * c[0][1] + c[0][2] * c[0][2] + c[1][2] * c[1][2];
   const double q = tr / 3;
   const double p2 = (c[0][0] - q) * (c[0][0] - q) + (c[1][1] - q) * (c[1][1] - q) + (c[2][2] - q) * (c[2][2] - q) + 2 * p1;
   const double p = std::sqrt(p2 / 6);
   if (p < 1e-12)
      return 1.0 / 3; // isotropic
   double bm[3][3];
   for (int a = 0; a < 3; ++a)
      for (int b = 0; b < 3; ++b)
         bm[a][b] = (c[a][b] - (a == b ? q : 0)) / p;
   const double detB = bm[0][0] * (bm[1][1] * bm[2][2] - bm[1][2] * bm[2][1]) -
                       bm[0][1] * (bm[1][0] * bm[2][2] - bm[1][2] * bm[2][0]) +
                       bm[0][2] * (bm[1][0] * bm[2][1] - bm[1][1] * bm[2][0]);
   const double r = std::max(-1.0, std::min(1.0, detB / 2));
   const double lMax = q + 2 * p * std::cos(std::acos(r) / 3);
   return lMax / tr;
}

/// hits in the DIGI frame; z_lab = zSign*z + zOff
inline Result estimate(std::vector<P> h, double B = 2.85, double zSign = -1, double zOff = 1000, int minPts = 5,
                       double binW = 4.0, double maxEdge = -4.0, int minBlob = 8, double blobFrac = 0.15,
                       double maxDropFrac = 0.5, EndOpt eo = EndOpt())
{
   Result r;
   for (size_t i = 0; i < h.size(); ++i) { h[i].z = zSign * h[i].z + zOff; h[i].id = i; }
   r.n = h.size();
   // BLOB REMOVAL: find the pieces the MST falls into when edges > maxEdge are cut, drop only SMALL
   // pieces (< max(minBlob, blobFrac*n) hits), then build the backbone on everything else WITHOUT a
   // cut. Keeping only the largest piece instead chopped real tracks (steep tracks have >10 mm 3D
   // gaps between pads): median 9-26 % of hits kept in data, and differently in sim (2026-10-30).
   {
      std::vector<int> comp;
      if (maxEdge != 0 && maxEdge < 1e8 && pieces(h, maxEdge, comp) > 1) {
         std::map<int, int> sz;
         for (int c : comp) ++sz[c];
         const int keepMin = std::max(minBlob, (int)(blobFrac * h.size()));
         // drop small pieces SMALLEST FIRST, never more than maxDropFrac of the hits in total: on a
         // dense multi-turn spiral the cut can split real track, and an uncapped rule then threw away
         // 1770 of 2098 hits (run_0016 ev 240).
         std::vector<std::pair<int, int>> small; // (size, piece)
         std::map<int, std::vector<int>> members;
         if (eo.blobLinMin > 0)
            for (size_t i = 0; i < h.size(); ++i) members[comp[i]].push_back(i);
         // circle of the LARGEST piece, for the on-circle test
         double cxL = 0, cyL = 0, RL = 0;
         bool haveL = false;
         if (eo.blobOffCircle > 0) {
            int cL = -1, nL = 0;
            for (auto &[c, n] : sz) if (n > nL) { nL = n; cL = c; }
            std::vector<int> idxL;
            for (size_t i = 0; i < h.size(); ++i) if (comp[i] == cL) idxL.push_back(i);
            std::vector<int> inlL;
            double cutL;
            haveL = robustCircle(h, idxL, cxL, cyL, RL, inlL, cutL);
         }
         for (auto &[c, n] : sz) {
            if (n >= keepMin) continue;
            if (haveL) { // on the track's circle -> chopped-off track, keep
               double dsum = 0;
               int nn = 0;
               for (size_t i = 0; i < h.size(); ++i)
                  if (comp[i] == c) { dsum += std::abs(std::hypot(h[i].x - cxL, h[i].y - cyL) - RL); ++nn; }
               if (nn && dsum / nn < eo.blobOffCircle) continue;
            }
            // a LINE-LIKE small piece is chopped-off track (steep, sparse forward protons), not a blob
            if (eo.blobLinMin > 0 && n >= 3 && linearity(h, members[c]) > eo.blobLinMin) continue;
            small.push_back({n, c});
         }
         std::sort(small.begin(), small.end());
         std::set<int> drop;
         int nDrop = 0;
         for (auto &[n, c] : small) {
            if (nDrop + n > maxDropFrac * h.size()) break;
            drop.insert(c);
            nDrop += n;
         }
         std::vector<P> kept;
         for (size_t i = 0; i < h.size(); ++i)
            if (!drop.count(comp[i])) kept.push_back(h[i]);
         if (kept.size() >= 5) { r.nDropped = h.size() - kept.size(); h.swap(kept); }
      }
   }
   if (!backbone(h, r.order, r.nBackbone, r.length, &r.s)) return r;
   r.hits = h;
   // CENTRELINE: charge-weighted centroids of the hits in binW-mm bins of backbone length. A track is
   // 3-4 pads WIDE in a1975 data; fitting the raw band biases short arcs to small circles.
   std::vector<P> c;
   std::vector<double> sc;
   for (size_t i = 0; i < r.order.size();) {
      size_t j = i;
      double wx = 0, wy = 0, wz = 0, q = 0, ws = 0, w = 0;
      while (j < r.order.size() && r.s[j] - r.s[i] < binW) {
         auto &p = h[r.order[j]];
         double wt = p.q > 0 ? p.q : 1;
         wx += wt * p.x; wy += wt * p.y; wz += wt * p.z; ws += wt * r.s[j]; w += wt; q += p.q;
         ++j;
      }
      c.push_back({wx / w, wy / w, wz / w, q});
      sc.push_back(ws / w);
      i = j;
   }
   r.cen = c;
   r.sc = sc;
   std::vector<int> oc(c.size());
   for (size_t i = 0; i < c.size(); ++i) oc[i] = i;
   r.e[0] = endHelix(c, oc, sc, B, minPts, eo.minLen, eo.rExcl, eo.guardFrac, eo.dcaMax, eo.rStartAxis);
   std::vector<int> rev(oc.rbegin(), oc.rend());
   std::vector<double> sr(sc.size());
   for (size_t i = 0; i < sc.size(); ++i) sr[i] = r.length - sc[sc.size() - 1 - i];
   r.e[1] = endHelix(c, rev, sr, B, minPts, eo.minLen, eo.rExcl, eo.guardFrac, eo.dcaMax, eo.rStartAxis);
   r.ok = r.e[0].ok && r.e[1].ok;
   // order -> indices of the CALLER's array (blob removal re-indexed the working copy)
   for (auto &o : r.order) o = h[o].id;
   return r;
}

// ---------------------------------------------------------------------------------------------
// VERTEX END by 2-of-3 majority (2026-10-30). Three independent calls:
//   sense   rotation sense seen travelling inward from an end; the vertex end shows `expSense`
//           (+1 measured in sim AND a1975 data, raw x-y). Abstains when both ends show the same sense.
//   spyral  caller-supplied end index from AtSpyralPID (or -1). Spyral forward = vertex at the LOWEST
//           z_DIGI end -- map it with the digi-frame z of the two ends, not the lab one.
//   bragg   charge per centreline bin on the small pads (rho <= rPad), last third / first third;
//           a rise marks the STOPPING end. Votes only if the track is contained (rho_max < rCont),
//           where it scored ~94 % on data; it has no peak to see on tracks that leave.
// 3/3 or 2/3 decide. One voter alone decides with nFor = 1. A 1-vs-1 split or no vote is AMBIGUOUS:
// the sense decides if it voted, else the helix closest approach, else end A.
struct EndVote {
   int end{0};
   int nFor{0}, nVotes{0};
   bool ambiguous{true};
   int sense{-1}, spyral{-1}, bragg{-1};
   double braggRatio{0};
};

inline int braggVote(const Result &r, double thr, double rPad, double rCont, double &ratio)
{
   ratio = 0;
   std::vector<double> qi;
   double rhoMax = 0;
   for (auto &c : r.cen) {
      double rho = std::hypot(c.x, c.y);
      rhoMax = std::max(rhoMax, rho);
      if (rho <= rPad) qi.push_back(c.q);
   }
   if (qi.size() < 6) return -1;
   size_t k = qi.size() / 3;
   double a = 0, b = 0;
   for (size_t i = 0; i < k; ++i) { a += qi[i]; b += qi[qi.size() - 1 - i]; }
   ratio = a > 0 ? b / a : 0;
   if (rhoMax >= rCont) return -1;
   if (ratio > 1 + thr) return 0;            // rises toward B: stops at B, vertex A
   if (ratio < 1 / (1 + thr)) return 1;      // rises toward A
   return -1;
}

inline EndVote majorityEnd(const Result &r, int spyralEnd, int expSense = 1, double braggThr = 0.2,
                           double rPad = 152.0, double rCont = 150.0)
{
   EndVote V;
   if (!r.ok) return V;
   V.sense = r.e[0].sense == r.e[1].sense ? -1 : (r.e[0].sense == expSense ? 0 : 1);
   V.spyral = spyralEnd;
   V.bragg = braggVote(r, braggThr, rPad, rCont, V.braggRatio);
   int n0 = 0, n1 = 0;
   for (int v : {V.sense, V.spyral, V.bragg}) {
      if (v == 0) ++n0;
      if (v == 1) ++n1;
   }
   V.nVotes = n0 + n1;
   if (n0 != n1) {
      V.end = n0 > n1 ? 0 : 1;
      V.nFor = std::max(n0, n1);
      V.ambiguous = false;
      return V;
   }
   V.ambiguous = true; // tie (1-1) or no vote
   if (V.sense >= 0) V.end = V.sense;
   else V.end = r.e[0].dca <= r.e[1].dca ? 0 : 1;
   V.nFor = V.nVotes ? 1 : 0;
   return V;
}

// ---------------------------------------------------------------------------------------------
// PID from one end (the estimator validated on a1975 D2, 2026-10-30): Brho and theta from that end's
// helix; sqrt(dE/dx) = summed centroid charge over the FIRST ARC (from the first inlier to the
// inlier farthest from it in x-y) restricted to the small pads (rho <= rPad, stopping at the first
// exit, as Spyral does), divided by the HELIX path of that segment sqrt(R^2+(dz/dphi)^2)*|dphi|.
// The centroid polyline length zigzags across a 3-4-pad-wide track and was the dominant dE/dx noise.
inline bool firstArcPID(const Result &r, int end, double rPad, double &brho, double &sq, double &theta,
                        double &path)
{
   const auto &E = r.e[end];
   if (!E.ok || E.brho <= 0)
      return false;
   const int m = r.cen.size();
   auto inl = [&](const P &p) { return std::abs(std::hypot(p.x - E.xc, p.y - E.yc) - E.R) < E.resCut; };
   auto at = [&](int jj) { return end == 0 ? jj : m - 1 - jj; };
   int j0 = at(0), jFar = 0;
   for (int jj = 0; jj < m; ++jj)
      if (inl(r.cen[at(jj)])) {
         j0 = at(jj);
         break;
      }
   double dFar = -1;
   for (int jj = 0; jj < m; ++jj) {
      const auto &p = r.cen[at(jj)];
      if (!inl(p))
         continue;
      double d = std::hypot(p.x - r.cen[j0].x, p.y - r.cen[j0].y);
      if (d > dFar) {
         dFar = d;
         jFar = jj;
      }
   }
   double q = 0, phPrev = 0, phAcc = 0, phFirst = 0, phLast = 0;
   bool started = false;
   for (int jj = 0; jj <= jFar; ++jj) {
      const auto &p = r.cen[at(jj)];
      double a = std::atan2(p.y - E.yc, p.x - E.xc);
      if (jj == 0)
         phPrev = a;
      else {
         double d = a - phPrev;
         while (d > M_PI)
            d -= 2 * M_PI;
         while (d < -M_PI)
            d += 2 * M_PI;
         phAcc += d;
         phPrev = a;
      }
      bool inner = std::hypot(p.x, p.y) <= rPad;
      if (!inl(p))
         continue;
      if (!inner && started)
         break;
      if (inner) {
         if (!started) {
            phFirst = phAcc;
            started = true;
         }
         phLast = phAcc;
         q += p.q;
      }
   }
   path = std::sqrt(E.R * E.R + E.dzdphi * E.dzdphi) * std::abs(phLast - phFirst);
   if (!started || path < 5)
      return false;
   brho = E.brho;
   sq = std::sqrt(q / path);
   theta = E.theta;
   return true;
}


/// TRUNCATED-MEAN variant of firstArcPID: the same first-arc, small-pad, inlier segment of the vertex-end
/// helix is cut into consecutive samples of at least sampleLen mm of HELIX path; each sample's dE/dx =
/// its charge / its helix path; sqrt(dE/dx) is taken from the mean of the lowest keepFrac of the
/// samples, which removes the Landau tail (the classic TPC estimator). keepFrac >= 1 reproduces the
/// plain charge-over-path ratio only approximately (samples are averaged, not pooled).
inline bool firstArcPIDTrunc(const Result &r, int end, double rPad, double sampleLen, double keepFrac, double &brho,
                             double &sq, double &theta, double &path, int &nSamples)
{
   const auto &E = r.e[end];
   nSamples = 0;
   if (!E.ok || E.brho <= 0)
      return false;
   const int m = r.cen.size();
   auto inl = [&](const P &p) { return std::abs(std::hypot(p.x - E.xc, p.y - E.yc) - E.R) < E.resCut; };
   auto at = [&](int jj) { return end == 0 ? jj : m - 1 - jj; };
   int j0 = at(0), jFar = 0;
   for (int jj = 0; jj < m; ++jj)
      if (inl(r.cen[at(jj)])) {
         j0 = at(jj);
         break;
      }
   double dFar = -1;
   for (int jj = 0; jj < m; ++jj) {
      const auto &p = r.cen[at(jj)];
      if (!inl(p))
         continue;
      double d = std::hypot(p.x - r.cen[j0].x, p.y - r.cen[j0].y);
      if (d > dFar) {
         dFar = d;
         jFar = jj;
      }
   }
   const double pathPerRad = std::sqrt(E.R * E.R + E.dzdphi * E.dzdphi);
   std::vector<double> samp;
   double qS = 0, phS = 0, phPrev = 0, phAcc = 0, total = 0;
   bool started = false;
   for (int jj = 0; jj <= jFar; ++jj) {
      const auto &p = r.cen[at(jj)];
      double a = std::atan2(p.y - E.yc, p.x - E.xc);
      if (jj == 0)
         phPrev = a;
      else {
         double d = a - phPrev;
         while (d > M_PI)
            d -= 2 * M_PI;
         while (d < -M_PI)
            d += 2 * M_PI;
         phAcc += d;
         phPrev = a;
      }
      bool inner = std::hypot(p.x, p.y) <= rPad;
      if (!inl(p))
         continue;
      if (!inner && started)
         break;
      if (!inner)
         continue;
      if (!started) {
         started = true;
         phS = phAcc;
      }
      qS += p.q;
      double len = pathPerRad * std::abs(phAcc - phS);
      if (len >= sampleLen) {
         samp.push_back(qS / len);
         total += len;
         qS = 0;
         phS = phAcc;
      }
   }
   nSamples = samp.size();
   if (nSamples < 3)
      return false;
   std::sort(samp.begin(), samp.end());
   size_t k = std::max<size_t>(1, (size_t)std::lround(keepFrac * samp.size()));
   k = std::min(k, samp.size());
   double sum = 0;
   for (size_t i = 0; i < k; ++i)
      sum += samp[i];
   brho = E.brho;
   sq = std::sqrt(sum / k);
   theta = E.theta;
   path = total;
   return true;
}
} // namespace Helix3D
} // namespace AtTools
#endif
