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

// Faithful HDBSCAN*: returns a label per point (-1 = noise).
std::vector<int> runHDBSCAN(const std::vector<std::array<double, 3>> &P, int minClusterSize, int minSamples,
                            bool allowSingleCluster)
{
   const int n = static_cast<int>(P.size());
   std::vector<int> labels(n, -1);
   if (n < minClusterSize || n < 2)
      return labels;

   auto dist = [&](int i, int j) {
      double dx = P[i][0] - P[j][0], dy = P[i][1] - P[j][1], dz = P[i][2] - P[j][2];
      return std::sqrt(dx * dx + dy * dy + dz * dz);
   };

   // 1) core distance = (minSamples)-th nearest-neighbour distance (incl. self at 0)
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

   // 2) MST of the mutual-reachability graph (dense Prim)
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

   // 3) single-linkage dendrogram from sorted MST (nodes 0..n-1 leaves, n..2n-2 internal)
   std::sort(mst.begin(), mst.end(), [](const E &x, const E &y) { return x.w < y.w; });
   const int nNodes = 2 * n - 1;
   std::vector<int> childL(nNodes, -1), childR(nNodes, -1);
   std::vector<int> nodeSize(nNodes, 1);
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
         int ra = froot(e.a), rb = froot(e.b);
         int na = node[ra], nb = node[rb];
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

   // leaves under a node (memoised lazily)
   std::function<void(int, std::vector<int> &)> gatherLeaves = [&](int node, std::vector<int> &out) {
      if (node < n) {
         out.push_back(node);
         return;
      }
      gatherLeaves(childL[node], out);
      gatherLeaves(childR[node], out);
   };

   // 4) condense the tree (min_cluster_size). Cluster ids start at n.
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
      int node = bfs[qi];
      if (node < n)
         continue;
      double lambda = nodeDist[node] > 0 ? 1.0 / nodeDist[node] : kInf;
      int L = childL[node], R = childR[node];
      int sL = nodeSize[L], sR = nodeSize[R];
      bool bigL = sL >= minClusterSize, bigR = sR >= minClusterSize;
      if (bigL && bigR) { // genuine split into two clusters
         relabel[L] = nextLabel++;
         cond.push_back({relabel[node], relabel[L], lambda, sL});
         bfs.push_back(L);
         relabel[R] = nextLabel++;
         cond.push_back({relabel[node], relabel[R], lambda, sR});
         bfs.push_back(R);
      } else if (!bigL && !bigR) { // cluster dies -> all its points fall to noise
         std::vector<int> lv;
         gatherLeaves(node, lv);
         for (int p : lv)
            cond.push_back({relabel[node], p, lambda, 1});
      } else { // one child keeps the cluster, the small one's points fall out
         int keep = bigL ? L : R, drop = bigL ? R : L;
         relabel[keep] = relabel[node];
         bfs.push_back(keep);
         std::vector<int> lv;
         gatherLeaves(drop, lv);
         for (int p : lv)
            cond.push_back({relabel[node], p, lambda, 1});
      }
   }

   // 5) cluster stability = sum over child rows of size*(lambda - lambda_birth)
   std::map<int, double> birth, stability;
   birth[n] = 0.0;
   for (const auto &r : cond)
      if (r.child >= n)
         birth[r.child] = r.lambda;
   for (const auto &r : cond)
      stability[r.parent] += r.size * (r.lambda - birth[r.parent]);

   // 6) excess-of-mass selection (process clusters bottom-up = descending id)
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
      isCluster[n] = false; // forbid the all-points root cluster
   for (int c : clusterIds) {
      if (!allowSingleCluster && c == n)
         continue;
      double childSum = 0;
      for (int cc : childClusters[c])
         childSum += stability[cc];
      if (childSum > stability[c]) {
         isCluster[c] = false;
         stability[c] = childSum; // propagate the better children value upward
      } else {
         // select c, unselect every descendant cluster
         std::vector<int> dq{c};
         for (size_t i = 0; i < dq.size(); ++i)
            for (int cc : childClusters[dq[i]]) {
               isCluster[cc] = false;
               dq.push_back(cc);
            }
      }
   }

   // 7) assign labels: every point under a selected cluster gets that cluster's label
   std::map<int, std::vector<int>> rowsByParent;
   for (size_t i = 0; i < cond.size(); ++i)
      rowsByParent[cond[i].parent].push_back(static_cast<int>(i));
   int lab = 0;
   for (int c : clusterIds) {
      if (!isCluster[c])
         continue;
      std::vector<int> stack{c};
      while (!stack.empty()) {
         int cur = stack.back();
         stack.pop_back();
         for (int ri : rowsByParent[cur]) {
            const auto &r = cond[ri];
            if (r.child < n)
               labels[r.child] = lab;
            else
               stack.push_back(r.child);
         }
      }
      ++lab;
   }
   return labels;
}

} // namespace

std::unique_ptr<AtPatternEvent> AtPATTERN::AtTrackFinderHDBSCAN::FindTracks(AtEvent &event)
{
   const int nHits = event.GetNumHits();
   std::vector<std::array<double, 3>> pts;
   pts.reserve(nHits);
   for (int i = 0; i < nHits; ++i) {
      auto p = event.GetHit(i).GetPosition();
      pts.push_back({p.X(), p.Y(), p.Z()});
   }

   std::vector<int> labels = runHDBSCAN(pts, fMinClusterSize, fMinSamples, fAllowSingleCluster);

   // group point indices by cluster label
   std::map<int, std::vector<int>> clusters;
   for (int i = 0; i < nHits; ++i)
      if (labels[i] >= 0)
         clusters[labels[i]].push_back(i);

   auto retEvent = std::make_unique<AtPatternEvent>();
   int trackID = 0;
   for (auto &kv : clusters) {
      AtTrack track;
      for (int pi : kv.second)
         track.AddHit(event.GetHit(pi));
      track.SetTrackID(trackID++);
      if (track.GetHitArray().size() > 0)
         SetTrackInitialParameters(track);
      retEvent->AddTrack(std::move(track));
   }
   for (int i = 0; i < nHits; ++i)
      if (labels[i] < 0)
         retEvent->AddNoise(event.GetHit(i));

   return retEvent;
}

ClassImp(AtPATTERN::AtTrackFinderHDBSCAN);
