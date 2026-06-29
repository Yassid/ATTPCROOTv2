/********************************************************************
 * AtTrackFinderHDBSCAN                                              *
 * Density-based clustering (HDBSCAN*) + Spyral-style cluster joins. *
 * Separates OVERLAPPING/crossing tracks where triplet clustering    *
 * (AtTrackFinderTC) merges them, then JOINS over-segmented pieces   *
 * of the same track (circle-overlap and/or end-to-end continuity).  *
 * Self-contained C++: core distances -> mutual-reachability MST ->  *
 * condensed tree -> excess-of-mass (with cluster_selection_epsilon) *
 * -> circle-fit cluster joining. Mirrors Spyral a1975 config.       *
 ********************************************************************/
#ifndef ATTRACKFINDERHDBSCAN_H
#define ATTRACKFINDERHDBSCAN_H

#include "AtPRA.h" // for AtPRA

#include <Rtypes.h> // for ClassDefOverride

#include <memory> // for unique_ptr
#include <string>

class AtEvent;
class AtPatternEvent;
class TBuffer;
class TClass;
class TMemberInspector;

namespace AtPATTERN {

class AtTrackFinderHDBSCAN : public AtPRA {
private:
   // --- HDBSCAN* core (Spyral a1975: min_points=3, adaptive min_size, cse=10) ---
   int fMinClusterSize{0};        ///< fixed min cluster size; if <=0 use adaptive max(lowerCutoff, scale*n)
   int fMinSamples{3};            ///< k for core distance (Spyral min_points)
   double fMinSizeScaleFactor{0.05}; ///< adaptive min_cluster_size = max(lowerCutoff, scale*nPoints)
   int fMinSizeLowerCutoff{10};
   double fClusterSelectionEpsilon{10.0}; ///< don't split clusters separated by < this (mm). Big over-seg lever.
   bool fAllowSingleCluster{true};

   // --- post-HDBSCAN joining (Spyral a1975: overlap join, ratio 0.25, min_size_join 15) ---
   std::string fJoinMethod{"overlap"}; ///< "overlap" | "continuity" | "none"
   double fCircleOverlapRatio{0.25};   ///< merge if circle area_overlap > ratio * smaller_area
   int fMinClusterSizeJoin{15};
   double fJoinZFraction{0.2};       ///< continuity: z-gap < join_z_fraction * (total z extent)
   double fJoinRadiusFraction{0.3};  ///< continuity: rho diff < join_radius_fraction * avg radius
   double fDirectionThreshold{0.5};  ///< forward/backward decision; only same-direction clusters merge
   // --- "motion" join: vertex-seeded smooth-trajectory following (lightweight equation-of-motion) ---
   double fMotionGapTol{40.0};       ///< max position gap (mm) between chunk ends to chain
   double fMotionAngleTol{35.0};     ///< max tangent-continuity angle (deg) at the junction
   // --- per-cluster Local Outlier Factor cleaning (Spyral outlier_scale_factor), after the join ---
   double fOutlierScaleFactor{0.05}; ///< LOF n_neighbors = scale*clusterSize (min 2); <=0 disables
   double fOutlierThreshold{1.5};    ///< LOF > this => outlier -> noise (sklearn contamination='auto')

public:
   AtTrackFinderHDBSCAN() = default;
   ~AtTrackFinderHDBSCAN() = default;

   std::unique_ptr<AtPatternEvent> FindTracks(AtEvent &event) override;

   void SetMinClusterSize(int n) { fMinClusterSize = n; }
   void SetMinSamples(int n) { fMinSamples = n; }
   void SetMinSizeScaleFactor(double f) { fMinSizeScaleFactor = f; }
   void SetMinSizeLowerCutoff(int n) { fMinSizeLowerCutoff = n; }
   void SetClusterSelectionEpsilon(double e) { fClusterSelectionEpsilon = e; }
   void SetAllowSingleCluster(bool v) { fAllowSingleCluster = v; }
   void SetJoinMethod(const std::string &m) { fJoinMethod = m; }
   void SetCircleOverlapRatio(double r) { fCircleOverlapRatio = r; }
   void SetMinClusterSizeJoin(int n) { fMinClusterSizeJoin = n; }
   void SetJoinZFraction(double f) { fJoinZFraction = f; }
   void SetJoinRadiusFraction(double f) { fJoinRadiusFraction = f; }
   void SetDirectionThreshold(double t) { fDirectionThreshold = t; }
   void SetMotionGapTol(double g) { fMotionGapTol = g; }
   void SetMotionAngleTol(double a) { fMotionAngleTol = a; }
   void SetOutlierScaleFactor(double s) { fOutlierScaleFactor = s; } // <=0 disables LOF cleaning
   void SetOutlierThreshold(double t) { fOutlierThreshold = t; }

   ClassDefOverride(AtTrackFinderHDBSCAN, 1);
};

} // namespace AtPATTERN

#endif
