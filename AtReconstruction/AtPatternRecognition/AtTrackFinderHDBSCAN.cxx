#include "AtTrackFinderHDBSCAN.h"

#include "AtEvent.h"
#include "AtHit.h"
#include "AtPatternEvent.h"
#include "AtTrack.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <vector>

namespace {

constexpr double kInf = std::numeric_limits<double>::infinity();
using P3 = std::array<double, 3>;

// ============================ HDBSCAN* core ============================
// Returns a label per point (-1 = noise). cse = cluster_selection_epsilon (mm): clusters separated
// by less than cse are not split (merged up the condensed tree).
std::vector<int> runHDBSCAN(const std::vector<P3> &P, int minClusterSize, int minSamples, bool allowSingleCluster,
                            double cse)
{
   const int n = static_cast<int>(P.size());
   std::vector<int> labels(n, -1);
   if (n < minClusterSize || n < 2)
      return labels;

   auto dist = [&](int i, int j) {
      double dx = P[i][0] - P[j][0], dy = P[i][1] - P[j][1], dz = P[i][2] - P[j][2];
      return std::sqrt(dx * dx + dy * dy + dz * dz);
   };

   const int k = std::max(1, std::min(minSamples, n - 1));
   std::vector<double> core(n);
   {
      std::vector<double> d(n);
      for (int i = 0; i < n; ++i) {
         for (int j = 0; j < n; ++j)
            d[j] = dist(i, j);
         std::nth_element(d.begin(), d.begin() + k, d.end());
         core[i] = d[k];
      }
   }
   auto mrd = [&](int i, int j) { return std::max(std::max(core[i], core[j]), dist(i, j)); };

   struct E {
      int a, b;
      double w;
   };
   std::vector<E> mst;
   mst.reserve(n - 1);
   {
      std::vector<char> inMST(n, 0);
      std::vector<double> key(n, kInf);
      std::vector<int> par(n, -1);
      key[0] = 0;
      for (int it = 0; it < n; ++it) {
         int u = -1;
         double best = kInf;
         for (int v = 0; v < n; ++v)
            if (!inMST[v] && key[v] < best) {
               best = key[v];
               u = v;
            }
         if (u < 0)
            break;
         inMST[u] = 1;
         if (par[u] >= 0)
            mst.push_back({par[u], u, key[u]});
         for (int v = 0; v < n; ++v)
            if (!inMST[v]) {
               double w = mrd(u, v);
               if (w < key[v]) {
                  key[v] = w;
                  par[v] = u;
               }
            }
      }
   }

   std::sort(mst.begin(), mst.end(), [](const E &x, const E &y) { return x.w < y.w; });
   const int nNodes = 2 * n - 1;
   std::vector<int> childL(nNodes, -1), childR(nNodes, -1), nodeSize(nNodes, 1);
   std::vector<double> nodeDist(nNodes, 0.0);
   {
      std::vector<int> rep(n);
      std::iota(rep.begin(), rep.end(), 0);
      std::function<int(int)> froot = [&](int x) {
         while (rep[x] != x) {
            rep[x] = rep[rep[x]];
            x = rep[x];
         }
         return x;
      };
      std::vector<int> node(n);
      std::iota(node.begin(), node.end(), 0);
      int nextId = n;
      for (const auto &e : mst) {
         int ra = froot(e.a), rb = froot(e.b), na = node[ra], nb = node[rb];
         childL[nextId] = na;
         childR[nextId] = nb;
         nodeDist[nextId] = e.w;
         nodeSize[nextId] = nodeSize[na] + nodeSize[nb];
         rep[ra] = rb;
         node[rb] = nextId;
         ++nextId;
      }
   }
   const int root = nNodes - 1;
   std::function<void(int, std::vector<int> &)> gatherLeaves = [&](int nd, std::vector<int> &out) {
      if (nd < n) { out.push_back(nd); return; }
      gatherLeaves(childL[nd], out);
      gatherLeaves(childR[nd], out);
   };

   struct CRow {
      int parent, child;
      double lambda;
      int size;
   };
   std::vector<CRow> cond;
   std::vector<int> relabel(nNodes, -1);
   relabel[root] = n;
   int nextLabel = n + 1;
   std::vector<int> bfs{root};
   for (size_t qi = 0; qi < bfs.size(); ++qi) {
      int nd = bfs[qi];
      if (nd < n)
         continue;
      double lambda = nodeDist[nd] > 0 ? 1.0 / nodeDist[nd] : kInf;
      int L = childL[nd], R = childR[nd], sL = nodeSize[L], sR = nodeSize[R];
      bool bigL = sL >= minClusterSize, bigR = sR >= minClusterSize;
      if (bigL && bigR) {
         relabel[L] = nextLabel++;
         cond.push_back({relabel[nd], relabel[L], lambda, sL});
         bfs.push_back(L);
         relabel[R] = nextLabel++;
         cond.push_back({relabel[nd], relabel[R], lambda, sR});
         bfs.push_back(R);
      } else if (!bigL && !bigR) {
         std::vector<int> lv;
         gatherLeaves(nd, lv);
         for (int p : lv)
            cond.push_back({relabel[nd], p, lambda, 1});
      } else {
         int keep = bigL ? L : R, drop = bigL ? R : L;
         relabel[keep] = relabel[nd];
         bfs.push_back(keep);
         std::vector<int> lv;
         gatherLeaves(drop, lv);
         for (int p : lv)
            cond.push_back({relabel[nd], p, lambda, 1});
      }
   }

   std::map<int, double> birth, stability;
   std::map<int, int> clusterParent;
   birth[n] = 0.0;
   for (const auto &r : cond)
      if (r.child >= n) {
         birth[r.child] = r.lambda;
         clusterParent[r.child] = r.parent;
      }
   for (const auto &r : cond)
      stability[r.parent] += r.size * (r.lambda - birth[r.parent]);

   std::map<int, std::vector<int>> childClusters;
   for (const auto &r : cond)
      if (r.child >= n)
         childClusters[r.parent].push_back(r.child);
   std::vector<int> clusterIds;
   for (auto &kv : stability)
      clusterIds.push_back(kv.first);
   std::sort(clusterIds.begin(), clusterIds.end(), std::greater<int>());

   std::map<int, bool> isCluster;
   for (int c : clusterIds)
      isCluster[c] = true;
   if (!allowSingleCluster)
      isCluster[n] = false;
   for (int c : clusterIds) {
      if (!allowSingleCluster && c == n)
         continue;
      double childSum = 0;
      for (int cc : childClusters[c])
         childSum += stability[cc];
      if (childSum > stability[c]) {
         isCluster[c] = false;
         stability[c] = childSum;
      } else {
         std::vector<int> dq{c};
         for (size_t i = 0; i < dq.size(); ++i)
            for (int cc : childClusters[dq[i]]) {
               isCluster[cc] = false;
               dq.push_back(cc);
            }
      }
   }

   // cluster_selection_epsilon: walk each selected cluster up while its birth split-distance < cse
   auto epsRep = [&](int c) {
      int w = c;
      while (clusterParent.count(w)) {
         double bornDist = birth[w] > 0 ? 1.0 / birth[w] : kInf;
         if (bornDist < cse)
            w = clusterParent[w];
         else
            break;
      }
      return w;
   };

   std::map<int, std::vector<int>> rowsByParent;
   for (size_t i = 0; i < cond.size(); ++i)
      rowsByParent[cond[i].parent].push_back(static_cast<int>(i));
   std::map<int, int> repLabel;
   int lab = 0;
   for (int c : clusterIds) {
      if (!isCluster[c])
         continue;
      int rep = (cse > 0) ? epsRep(c) : c;
      if (!repLabel.count(rep))
         repLabel[rep] = lab++;
      int myLab = repLabel[rep];
      std::vector<int> stack{c};
      while (!stack.empty()) {
         int cur = stack.back();
         stack.pop_back();
         for (int ri : rowsByParent[cur]) {
            const auto &r = cond[ri];
            if (r.child < n)
               labels[r.child] = myLab;
            else
               stack.push_back(r.child);
         }
      }
   }
   return labels;
}

// ============================ circle fit / direction / join ============================
// Analytic least-squares (Kasa) circle fit. Returns false if singular.
bool fitCircle(const std::vector<int> &idx, const std::vector<P3> &P, double &cx, double &cy, double &r)
{
   int n = idx.size();
   if (n < 3)
      return false;
   double mx = 0, my = 0;
   for (int i : idx) { mx += P[i][0]; my += P[i][1]; }
   mx /= n; my /= n;
   double Suu = 0, Svv = 0, Suv = 0, Suuu = 0, Svvv = 0, Suuv = 0, Suvv = 0;
   for (int i : idx) {
      double u = P[i][0] - mx, v = P[i][1] - my;
      Suu += u * u; Svv += v * v; Suv += u * v;
      Suuu += u * u * u; Svvv += v * v * v; Suuv += u * u * v; Suvv += u * v * v;
   }
   double det = Suu * Svv - Suv * Suv;
   if (std::fabs(det) < 1e-9)
      return false;
   double b1 = 0.5 * (Suuu + Suvv), b2 = 0.5 * (Svvv + Suuv);
   double uc = (b1 * Svv - b2 * Suv) / det, vc = (Suu * b2 - Suv * b1) / det;
   cx = uc + mx; cy = vc + my;
   r = std::sqrt(uc * uc + vc * vc + (Suu + Svv) / n);
   return true;
}

// Cluster winding direction vs z: +1 forward, -1 backward, 0 none (Spyral get_direction: sign-fraction
// of angle steps around the circle center, ordered by z).
int clusterDirection(std::vector<int> idx, const std::vector<P3> &P, double thr)
{
   if (idx.size() < 6)
      return 0;
   std::sort(idx.begin(), idx.end(), [&](int a, int b) { return P[a][2] < P[b][2]; });
   double cx, cy, r;
   if (!fitCircle(idx, P, cx, cy, r))
      return 0;
   std::vector<double> ang;
   ang.reserve(idx.size());
   for (int i : idx) {
      double a = std::atan2(P[i][1] - cy, P[i][0] - cx);
      if (a < 0) a += 2 * M_PI;
      ang.push_back(a);
   }
   int pos = 0, tot = 0;
   for (size_t i = 1; i < ang.size(); ++i) {
      if (ang[i] - ang[i - 1] > 0) pos++;
      tot++;
   }
   if (tot == 0) return 0;
   double pf = (double)pos / tot;
   if (pf > thr) return 1;
   if (1 - pf > thr) return -1;
   return 0;
}

double circleOverlapArea(double x1, double y1, double r1, double x2, double y2, double r2)
{
   double d = std::sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2));
   if (d >= r1 + r2)
      return 0.0;
   if (d <= std::fabs(r1 - r2))
      return M_PI * std::min(r1, r2) * std::min(r1, r2);
   double t1 = (d * d + r1 * r1 - r2 * r2) / (2 * d * r1);
   double t2 = (d * d + r2 * r2 - r1 * r1) / (2 * d * r2);
   t1 = std::min(1.0, std::max(-1.0, t1));
   t2 = std::min(1.0, std::max(-1.0, t2));
   double t3 = (-d + r1 + r2) * (d + r1 - r2) * (d - r1 + r2) * (d + r1 + r2);
   if (t3 < 0) t3 = 0;
   return r1 * r1 * std::acos(t1) + r2 * r2 * std::acos(t2) - 0.5 * std::sqrt(t3);
}

// merge clusters in-place (union-find per pass) until stable
void joinClusters(std::vector<std::vector<int>> &clusters, const std::vector<P3> &P, const std::string &method,
                  double ratio, int minSizeJoin, double zFrac, double rhoFrac, double dirThr)
{
   if (method == "none")
      return;
   bool changed = true;
   while (changed && clusters.size() > 1) {
      changed = false;
      int m = clusters.size();
      std::vector<double> cx(m), cy(m), r(m);
      std::vector<int> dir(m);
      std::vector<char> ok(m, 0);
      for (int i = 0; i < m; ++i) {
         if ((int)clusters[i].size() < minSizeJoin) continue;
         if (!fitCircle(clusters[i], P, cx[i], cy[i], r[i])) continue;
         if (std::isnan(cx[i]) || r[i] < 10.0) continue;
         dir[i] = clusterDirection(clusters[i], P, dirThr);
         ok[i] = 1;
      }
      std::vector<int> uf(m);
      std::iota(uf.begin(), uf.end(), 0);
      std::function<int(int)> find = [&](int x) { while (uf[x] != x) { uf[x] = uf[uf[x]]; x = uf[x]; } return x; };
      for (int i = 0; i < m; ++i) {
         if (!ok[i]) continue;
         for (int j = i + 1; j < m; ++j) {
            if (!ok[j] || dir[i] != dir[j]) continue;
            bool merge = false;
            if (method == "overlap") {
               double ov = circleOverlapArea(cx[i], cy[i], r[i], cx[j], cy[j], r[j]);
               double smaller = M_PI * std::min(r[i] * r[i], r[j] * r[j]);
               merge = ov > ratio * smaller;
            } else { // continuity: end-to-end on the same circle
               auto endRho = [&](const std::vector<int> &cl, double ccx, double ccy, bool up, double &zc) {
                  std::vector<int> s = cl;
                  std::sort(s.begin(), s.end(), [&](int a, int b) { return P[a][2] < P[b][2]; });
                  int k0 = up ? 0 : (int)s.size() - 3;
                  double sx = 0, sy = 0, sz = 0;
                  int cnt = 0;
                  for (int t = std::max(0, k0); t < std::min((int)s.size(), k0 + 3); ++t) {
                     sx += P[s[t]][0] - ccx; sy += P[s[t]][1] - ccy; sz += P[s[t]][2]; cnt++;
                  }
                  zc = sz / cnt;
                  return std::sqrt((sx / cnt) * (sx / cnt) + (sy / cnt) * (sy / cnt));
               };
               double zi, zj;
               double rhoi = endRho(clusters[i], cx[i], cy[i], false, zi); // i downstream end
               double rhoj = endRho(clusters[j], cx[j], cy[j], true, zj);  // j upstream end
               double zlo = 1e18, zhi = -1e18;
               for (int p : clusters[i]) { zlo = std::min(zlo, P[p][2]); zhi = std::max(zhi, P[p][2]); }
               for (int p : clusters[j]) { zlo = std::min(zlo, P[p][2]); zhi = std::max(zhi, P[p][2]); }
               double zThresh = (zhi - zlo) * zFrac;
               double rhoThresh = 0.5 * (r[i] + r[j]) * rhoFrac;
               merge = (std::fabs(zi - zj) < zThresh) && (std::fabs(rhoi - rhoj) < rhoThresh);
            }
            if (merge) {
               int ri = find(i), rj = find(j);
               if (ri != rj) { uf[rj] = ri; changed = true; }
            }
         }
      }
      if (changed) {
         std::map<int, std::vector<int>> groups;
         for (int i = 0; i < m; ++i)
            for (int p : clusters[i])
               groups[find(i)].push_back(p);
         clusters.clear();
         for (auto &kv : groups)
            clusters.push_back(std::move(kv.second));
      }
   }
}

// "motion" join: vertex-seeded smooth-trajectory following (lightweight equation-of-motion proxy).
// Anchors at the vertex (low-transverse-radius region), then chains chunks whose ends connect with
// CONTINUOUS tangent (direction of travel), allowing curvature to evolve (energy loss) -- unlike the
// fixed-circle joins. gapTol = max end-to-end gap (mm); angleTolDeg = max tangent mismatch at junction.
void motionJoin(std::vector<std::vector<int>> &clusters, const std::vector<P3> &P, double gapTol,
                double angleTolDeg, int minSizeJoin)
{
   int m = clusters.size();
   if (m < 2)
      return;
   const double cosTol = std::cos(angleTolDeg * M_PI / 180.0);
   auto d3 = [](const P3 &a, const P3 &b) {
      double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
      return std::sqrt(dx * dx + dy * dy + dz * dz);
   };
   auto dot = [](const P3 &a, const P3 &b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; };

   // vertex proxy: centroid of the lowest-transverse-radius decile of all clustered hits
   std::vector<double> rxy;
   for (auto &cl : clusters)
      for (int p : cl)
         rxy.push_back(std::sqrt(P[p][0] * P[p][0] + P[p][1] * P[p][1]));
   if (rxy.empty())
      return;
   std::vector<double> srt = rxy;
   std::sort(srt.begin(), srt.end());
   double rcut = srt[srt.size() / 10];
   P3 V{0, 0, 0};
   int vcnt = 0;
   for (auto &cl : clusters)
      for (int p : cl)
         if (std::sqrt(P[p][0] * P[p][0] + P[p][1] * P[p][1]) <= rcut) {
            V[0] += P[p][0]; V[1] += P[p][1]; V[2] += P[p][2]; vcnt++;
         }
   if (vcnt > 0) { V[0] /= vcnt; V[1] /= vcnt; V[2] /= vcnt; }

   // per-chunk: near/far ends (near = closer to V) + travel tangent at each (oriented near->far)
   struct Ends { P3 nearP, farP, nearT, farT; bool ok; };
   std::vector<Ends> E(m);
   for (int c = 0; c < m; ++c) {
      E[c].ok = false;
      auto &cl = clusters[c];
      if ((int)cl.size() < minSizeJoin)
         continue;
      int e0 = cl[0], e1 = cl[0];
      double best = -1;
      for (size_t i = 0; i < cl.size(); ++i)
         for (size_t j = i + 1; j < cl.size(); ++j) {
            double dd = d3(P[cl[i]], P[cl[j]]);
            if (dd > best) { best = dd; e0 = cl[i]; e1 = cl[j]; }
         }
      P3 p0 = P[e0], p1 = P[e1];
      if (d3(p0, V) > d3(p1, V)) std::swap(p0, p1);
      E[c].nearP = p0; E[c].farP = p1;
      auto travelTangent = [&](P3 endP, bool outward) {
         std::vector<std::pair<double, int>> dn;
         for (int p : cl) dn.push_back({d3(P[p], endP), p});
         std::sort(dn.begin(), dn.end());
         int K = std::min((int)dn.size(), 8);
         P3 cen{0, 0, 0};
         for (int t = 0; t < K; ++t) { cen[0] += P[dn[t].second][0]; cen[1] += P[dn[t].second][1]; cen[2] += P[dn[t].second][2]; }
         cen[0] /= K; cen[1] /= K; cen[2] /= K;
         P3 dir = outward ? P3{endP[0] - cen[0], endP[1] - cen[1], endP[2] - cen[2]}
                          : P3{cen[0] - endP[0], cen[1] - endP[1], cen[2] - endP[2]};
         double nn = std::sqrt(dot(dir, dir));
         if (nn > 1e-9) { dir[0] /= nn; dir[1] /= nn; dir[2] /= nn; }
         return dir;
      };
      E[c].nearT = travelTangent(p0, false); // travel direction entering the near end
      E[c].farT = travelTangent(p1, true);   // travel direction leaving the far end
      E[c].ok = true;
   }

   // greedy: seed from the chunk nearest the vertex, chain outward by tangent continuity
   std::vector<int> order(m);
   std::iota(order.begin(), order.end(), 0);
   std::sort(order.begin(), order.end(), [&](int a, int b) {
      double da = E[a].ok ? d3(E[a].nearP, V) : 1e18, db = E[b].ok ? d3(E[b].nearP, V) : 1e18;
      return da < db;
   });
   std::vector<int> group(m);
   std::iota(group.begin(), group.end(), 0);
   std::vector<char> used(m, 0);
   for (int s : order) {
      if (!E[s].ok || used[s])
         continue;
      used[s] = 1;
      P3 leadP = E[s].farP, leadT = E[s].farT;
      while (true) {
         int bestC = -1;
         double bestGap = gapTol;
         for (int c = 0; c < m; ++c) {
            if (!E[c].ok || used[c]) continue;
            double gap = d3(leadP, E[c].nearP);
            if (gap >= bestGap) continue;
            if (dot(leadT, E[c].nearT) < cosTol) continue; // tangent continuity
            P3 g{E[c].nearP[0] - leadP[0], E[c].nearP[1] - leadP[1], E[c].nearP[2] - leadP[2]};
            double gn = std::sqrt(dot(g, g));
            if (gn > 1e-6) { g[0] /= gn; g[1] /= gn; g[2] /= gn; if (dot(leadT, g) < cosTol) continue; } // gap points along travel
            bestC = c; bestGap = gap;
         }
         if (bestC < 0) break;
         used[bestC] = 1;
         group[bestC] = s;
         leadP = E[bestC].farP;
         leadT = E[bestC].farT;
      }
   }
   std::map<int, std::vector<int>> merged;
   for (int c = 0; c < m; ++c)
      for (int p : clusters[c])
         merged[group[c]].push_back(p);
   clusters.clear();
   for (auto &kv : merged)
      clusters.push_back(std::move(kv.second));
}

} // namespace

std::unique_ptr<AtPatternEvent> AtPATTERN::AtTrackFinderHDBSCAN::FindTracks(AtEvent &event)
{
   const int nHits = event.GetNumHits();
   std::vector<P3> pts;
   pts.reserve(nHits);
   for (int i = 0; i < nHits; ++i) {
      auto p = event.GetHit(i).GetPosition();
      pts.push_back({p.X(), p.Y(), p.Z()});
   }

   // adaptive min_cluster_size (Spyral) if fMinClusterSize <= 0
   int mcs = fMinClusterSize > 0
                ? fMinClusterSize
                : std::max(fMinSizeLowerCutoff, (int)(fMinSizeScaleFactor * nHits));
   std::vector<int> labels = runHDBSCAN(pts, mcs, fMinSamples, fAllowSingleCluster, fClusterSelectionEpsilon);

   std::map<int, std::vector<int>> byLabel;
   for (int i = 0; i < nHits; ++i)
      if (labels[i] >= 0)
         byLabel[labels[i]].push_back(i);
   std::vector<std::vector<int>> clusters;
   for (auto &kv : byLabel)
      clusters.push_back(std::move(kv.second));

   // Spyral-style joining of over-segmented pieces. "both" = continuity (safe z-stitching) then overlap.
   // "motion" = vertex-seeded smooth-trajectory following (lightweight equation-of-motion).
   if (fJoinMethod == "motion") {
      motionJoin(clusters, pts, fMotionGapTol, fMotionAngleTol, fMinClusterSizeJoin);
   } else if (fJoinMethod == "mover") { // motion (end-to-end energy-loss) then overlap (spiral splits)
      motionJoin(clusters, pts, fMotionGapTol, fMotionAngleTol, fMinClusterSizeJoin);
      joinClusters(clusters, pts, "overlap", fCircleOverlapRatio, fMinClusterSizeJoin, fJoinZFraction,
                   fJoinRadiusFraction, fDirectionThreshold);
   } else if (fJoinMethod == "both") {
      joinClusters(clusters, pts, "continuity", fCircleOverlapRatio, fMinClusterSizeJoin, fJoinZFraction,
                   fJoinRadiusFraction, fDirectionThreshold);
      joinClusters(clusters, pts, "overlap", fCircleOverlapRatio, fMinClusterSizeJoin, fJoinZFraction,
                   fJoinRadiusFraction, fDirectionThreshold);
   } else {
      joinClusters(clusters, pts, fJoinMethod, fCircleOverlapRatio, fMinClusterSizeJoin, fJoinZFraction,
                   fJoinRadiusFraction, fDirectionThreshold);
   }

   auto retEvent = std::make_unique<AtPatternEvent>();
   std::vector<char> claimed(nHits, 0);
   int trackID = 0;
   for (auto &cl : clusters) {
      AtTrack track;
      for (int pi : cl) { track.AddHit(event.GetHit(pi)); claimed[pi] = 1; }
      track.SetTrackID(trackID++);
      if (track.GetHitArray().size() > 0)
         SetTrackInitialParameters(track);
      retEvent->AddTrack(std::move(track));
   }
   for (int i = 0; i < nHits; ++i)
      if (!claimed[i])
         retEvent->AddNoise(event.GetHit(i));

   return retEvent;
}

ClassImp(AtPATTERN::AtTrackFinderHDBSCAN);
