/// @file AtPatternBridgeTask.h
/// @brief Fills the geometric/cluster fields on AtPatternEvent tracks that
/// were produced by AtSampleConsensusTask (RANSAC), so the downstream UKF
/// fitter can consume them like AtPRAtask output.
///
/// AtSampleConsensus stores hits + a fitted AtPattern but leaves
/// AtTrack::fGeoRadius / fGeoTheta / fGeoCenter at their defaults (zero) and
/// does not populate fHitClusterArray. The UKF fitter rejects tracks with
/// GeoRadius <= 0 and requires clusters. This task bridges the gap:
///   - For each track whose pattern is AtPatternCircle2D, copy
///     center/radius into Geo* fields and derive GeoTheta from the cluster
///     z-progression after Smooth3D clustering.
///   - Run AtTrackTransformer::ClusterizeSmooth3D to fill fHitClusterArray.

#ifndef ATPATTERNBRIDGETASK_H
#define ATPATTERNBRIDGETASK_H

#include <FairTask.h>

#include <Rtypes.h>
#include <TString.h>

class TBuffer;
class TClass;
class TMemberInspector;
class TClonesArray;

class AtPatternBridgeTask : public FairTask {
public:
   AtPatternBridgeTask() = default;
   ~AtPatternBridgeTask() override = default;

   InitStatus Init() override;
   void Exec(Option_t *option) override;

   void SetInputBranch(TString name) { fInputBranch = std::move(name); }
   void SetClusterRadius(double r_mm) { fClusterRadius = r_mm; }
   void SetClusterDistance(double d_mm) { fClusterDistance = d_mm; }

private:
   TString fInputBranch{"AtPatternEvent"};
   TClonesArray *fPatternEventArray{nullptr};
   double fClusterRadius{8.0};   // mm
   double fClusterDistance{5.0}; // mm

   ClassDefOverride(AtPatternBridgeTask, 1);
};

#endif
