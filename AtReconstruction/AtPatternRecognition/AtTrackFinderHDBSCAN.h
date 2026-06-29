/********************************************************************
 * AtTrackFinderHDBSCAN                                              *
 * Density-based clustering (HDBSCAN*) pattern recognition.          *
 * Separates OVERLAPPING/crossing tracks where triplet clustering    *
 * (AtTrackFinderTC) merges them. Self-contained C++ implementation: *
 * core distances -> mutual-reachability MST -> condensed tree ->    *
 * excess-of-mass (EOM) cluster extraction.                          *
 ********************************************************************/
#ifndef ATTRACKFINDERHDBSCAN_H
#define ATTRACKFINDERHDBSCAN_H

#include "AtPRA.h" // for AtPRA

#include <Rtypes.h> // for ClassDefOverride

#include <memory> // for unique_ptr

class AtEvent;
class AtPatternEvent;
class TBuffer;
class TClass;
class TMemberInspector;

namespace AtPATTERN {

class AtTrackFinderHDBSCAN : public AtPRA {
private:
   int fMinClusterSize{15};       ///< smallest accepted cluster (HDBSCAN min_cluster_size)
   int fMinSamples{6};            ///< k for core distance / how conservative the noise (HDBSCAN min_samples)
   bool fAllowSingleCluster{false}; ///< sklearn default: forbid the all-points root cluster -> best overlap
                                    ///< separation. Set true to let a lone track survive as one cluster.

public:
   AtTrackFinderHDBSCAN() = default;
   ~AtTrackFinderHDBSCAN() = default;

   std::unique_ptr<AtPatternEvent> FindTracks(AtEvent &event) override;

   void SetMinClusterSize(int n) { fMinClusterSize = n; }
   void SetMinSamples(int n) { fMinSamples = n; }
   void SetAllowSingleCluster(bool v) { fAllowSingleCluster = v; }

   ClassDefOverride(AtTrackFinderHDBSCAN, 1);
};

} // namespace AtPATTERN

#endif
