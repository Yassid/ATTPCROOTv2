#ifndef ATTRACKTRANSFORMER_H
#define ATTRACKTRANSFORMER_H

#include <Rtypes.h>

class AtTrack;

namespace AtTools {

class AtTrackTransformer {

public:
   AtTrackTransformer();
   ~AtTrackTransformer();

   /// Set diffusion and drift parameters for cluster covariance calculation.
   /// If not called, defaults are used (which may not match the detector).
   /// @param coefT Transverse diffusion coefficient [cm^2/us]
   /// @param coefL Longitudinal diffusion coefficient [cm^2/us]
   /// @param driftVel Electron drift velocity [cm/us]
   /// @param tbTime Time bucket duration [us]
   /// @param padResXY Pad position resolution in X/Y [mm] (default: padSize/sqrt(12))
   void SetDiffusionParams(double coefT, double coefL, double driftVel, double tbTime, double padResXY = -1);

   /// Per-hit transverse (xy) variance σ_xy²(z) in mm² for a hit at drift
   /// distance z_mm: padResXY² + 20·coefT·z_mm/vDrift. Matches the model
   /// used in BuildClusterFromHits and AtClusterize::getTransverseDiffusion.
   double GetSigmaXY2(double z_mm) const
   {
      double varT = (z_mm > 0) ? 20.0 * fCoefT * z_mm / fDriftVel : 0.0;
      return fPadResXY * fPadResXY + varT;
   }

   void ClusterizeSmooth3D(AtTrack &track, Float_t radius, Float_t distance);

   /// @brief Group consecutive hits into clusters of fixed size along the track.
   ///
   /// Simpler alternative to ClusterizeSmooth3D: no smoothing pass, no midpoint
   /// insertion. Groups `hitsPerCluster` consecutive hits and computes charge-weighted
   /// centroid. Produces clusters that sit directly on the track.
   /// @param track Track with raw hits
   /// @param hitsPerCluster Number of hits per cluster (default 15 → ~18 clusters for 280 hits)
   void ClusterizeByGroup(AtTrack &track, int hitsPerCluster = 15);

   /// @brief Arc-walk clustering: orders hits along a spanning tree arc, then groups
   /// consecutive hits into clusters. Gap-immune alternative to ClusterizeSmooth3D.
   ///
   /// Algorithm: build kNN graph from raw hits, find maximum spanning tree (Prim's),
   /// extract arc ordering via double-DFS, group hits into clusters with charge-weighted
   /// centroids and diffusion-based covariances.
   ///
   /// The number of hits per cluster is adaptive: nHits / targetClusters, with a
   /// minimum of minHitsPerCluster to ensure meaningful covariance estimation.
   /// @param track Track with raw hits (clusters written to track's cluster array)
   /// @param targetClusters Desired number of clusters per track (default 25)
   /// @param minHitsPerCluster Minimum hits per cluster floor (default 3)
   /// @param kNN Number of nearest neighbors for adjacency graph (default 10)
   void ClusterizeArcWalk(AtTrack &track, int targetClusters = 25, int minHitsPerCluster = 3, int kNN = 10);

   const std::tuple<Double_t, Double_t> GetPIDFromHits(AtTrack &track, Double_t theta);

   Bool_t FindVertexTrack(AtTrack *trA, AtTrack *trB);

   Bool_t MergeTracks(std::vector<AtTrack *> *trackCandSource, std::vector<AtTrack> *trackDest,
                      Bool_t enableSingleVertexTrack, Double_t clusterRadius, Double_t clusterDistance);

private:
   double fCoefT{0.00009};      ///< Transverse diffusion coefficient [cm^2/us]
   double fCoefL{0.0000009};    ///< Longitudinal diffusion coefficient [cm^2/us]
   double fDriftVel{1.0};       ///< Electron drift velocity [cm/us]
   double fTBTime{0.320};       ///< Time bucket duration [us]
   double fPadResXY{2.3};       ///< Pad position resolution in X/Y [mm] (default: 8mm/sqrt(12))
};

} // namespace AtTools

#endif
