#ifndef AtPRATASK_H
#define AtPRATASK_H

#include <FairTask.h> // for FairTask, InitStatus

#include <Rtypes.h>       // for Int_t, Double_t, Bool_t, THashConsistencyH...
#include <TClonesArray.h> // for TClonesArray
#include <TString.h>

#include <cstddef> // for size_t
#include <memory>  // for unique_ptr
#include <utility>

class AtDigiPar;
class TBuffer;
class TClass;
class TMemberInspector;
namespace AtPATTERN {
class AtPRA;
}

/**
 * @brief Task for finding patterns in hit clouds.
 *
 * Logic is in class AtPRA and derived types.
 */
class AtPRAtask : public FairTask {
private:
   TString fInputBranchName;
   TString fOutputBranchName;

   TClonesArray *fEventHArray{};
   TClonesArray fPatternEventArray;

   AtDigiPar *fPar;

   std::unique_ptr<AtPATTERN::AtPRA> fPRA; //<! Pattern-recognition algorithm (owned)
   bool fInjected{false};                  //<! true when the algorithm was supplied via the ctor

   Int_t fPRAlgorithm;

   Int_t fMinNumHits;
   Int_t fMaxNumHits;

   Bool_t kIsPersistence;

   // HC parameters
   float fHCs;
   size_t fHCk;
   size_t fHCn;
   size_t fHCm;
   float fHCr;
   float fHCa;
   float fHCt;
   size_t fHCpadding;

   // Prunning parameters
   Bool_t kSetPrunning;
   Int_t fKNN;             //<! Number of nearest neighbors kNN
   Double_t fStdDevMulkNN; //<! Std dev multiplier for kNN
   Double_t fkNNDist;      //<! Distance threshold for outlier rejection in kNN

   // Clustering parameters
   Double_t fClusterRadius{20.0};  // Overlapping clusters: radius > distance
   Double_t fClusterDistance{15.0}; // Gives 2.2% RMS (was r10 d20 → 3.8%)

   // Arc-walk clustering parameters
   bool fUseArcWalk{false};
   int fTargetClusters{25};
   int fMinHitsPerCluster{3};
   int fArcWalkKNN{10};

   // Drift-aware weighting for the circle radius fit (SetTrackInitialParameters)
   bool fUseDriftAwareWeights{false};
   double fPadResXYOverride{-1}; // mm, <0 keeps the auto-default
   bool fPreclusterRadiusFit{false};
   double fPreclusterBin_mm{3.0};

   // TC (fPRAlgorithm == 0) — primary/fragment heuristic switch (disable for annular geometries)
   bool fTCUseSelectAndMerge{true};
   bool fTCUseCircleMerge{false};
   double fTCCircleMergeRadiusTol{0.2};
   double fTCCircleMergeCenterTol{25.0};
   bool fChargeFromCenter{false};

   // Riemann parameters (fPRAlgorithm == 3)
   double fRiemannInlierDist{15.0};
   int fRiemannMinHits{10};
   int fRiemannMaxTracks{6};
   int fRiemannMaxIter{400};
   bool fRiemannUseSelectAndMerge{true};
   bool fRiemannUseArcWalkExtend{false};
   int fRiemannArcWalkWindow{10};
   int fRiemannArcWalkMaxMiss{5};

public:
   AtPRAtask();
   /// Dependency-injection constructor (mirrors AtFitterTask): drive a pre-configured
   /// pattern-recognition algorithm supplied by the caller. Configure the concrete
   /// AtPRA (e.g. AtTrackFinderTC / AtTrackFinderRiemann) fully before injecting; the
   /// task only forwards runtime AtDigiPar diffusion parameters and handles branch I/O.
   /// The legacy default constructor + SetPRAlgorithm() path is kept for existing macros.
   explicit AtPRAtask(std::unique_ptr<AtPATTERN::AtPRA> pra);
   ~AtPRAtask();

   virtual InitStatus Init();
   virtual void Exec(Option_t *option);
   virtual void SetParContainers();
   virtual void Finish();

   void SetInputBranch(TString branchName) { fInputBranchName = std::move(branchName); }
   void SetOutputBranch(TString branchName) { fOutputBranchName = std::move(branchName); }

   void SetPersistence(Bool_t value = kTRUE);
   void SetPRAlgorithm(Int_t value = 0);

   void SetScluster(float s) { fHCs = s; }
   void SetKtriplet(size_t k) { fHCk = k; }
   void SetNtriplet(size_t n) { fHCn = n; }
   void SetMcluster(size_t m) { fHCm = m; }
   void SetRsmooth(float r) { fHCr = r; }
   void SetAtriplet(float a) { fHCa = a; }
   void SetTcluster(float t) { fHCt = t; }
   void SetPadding(size_t padding) { fHCpadding = padding; }

   void SetMaxNumHits(Int_t maxHits) { fMaxNumHits = maxHits; }
   void SetMinNumHits(Int_t minHits) { fMinNumHits = minHits; }

   void SetPrunning() { kSetPrunning = kTRUE; }
   void SetkNN(Double_t knn) { fKNN = knn; }
   void SetStdDevMulkNN(Double_t stdDevMul) { fStdDevMulkNN = stdDevMul; }
   void SetkNNDist(Double_t dist) { fkNNDist = dist; }

   void SetClusterRadius(Double_t clusterRadius) { fClusterRadius = clusterRadius; }
   void SetClusterDistance(Double_t clusterDistance) { fClusterDistance = clusterDistance; }

   /// Use arc-walk clustering instead of ClusterizeSmooth3D (gap-immune)
   void SetUseArcWalk(bool use = true) { fUseArcWalk = use; }
   /// Target number of clusters per track (hits/cluster adapts automatically)
   void SetTargetClusters(int n) { fTargetClusters = n; }
   /// Minimum hits per cluster (floor for adaptive sizing)
   void SetMinHitsPerCluster(int n) { fMinHitsPerCluster = n; }
   /// Number of nearest neighbors for arc-walk kNN graph
   void SetArcWalkKNN(int k) { fArcWalkKNN = k; }

   /// Weight the circle radius fit by 1/σ²_xy(z) using diffusion params
   /// from AtDigiPar (padResXY² + 20·coefT·z/vDrift). Default false.
   void SetUseDriftAwareWeights(bool use = true) { fUseDriftAwareWeights = use; }
   /// Override the transverse pad position resolution (mm). Needed when the
   /// pad geometry differs from the AT-TPC default (e.g. pi-TPC's 2 mm
   /// square pads → padResXY = 2/√12 ≈ 0.58 mm; PUMA's annular ring varies).
   /// Negative value (default) keeps the built-in AtTrackTransformer default
   /// of 2.3 mm.
   void SetPadResXY(double mm) { fPadResXYOverride = mm; }
   /// Pre-cluster hits into a coarse 3D grid (charge-weighted centroids)
   /// before the PRA circle fit. Mitigates per-pad ±pad-pitch jitter from
   /// AtPSAMax in geometries with closely-spaced firing pads. Default false.
   void SetPreclusterRadiusFit(bool use = true) { fPreclusterRadiusFit = use; }
   void SetPreclusterBin(double mm) { fPreclusterBin_mm = mm; }

   /// TC (PRAlgorithm == 0): apply the AT-TPC primary/fragment heuristic. Disable for PUMA / annular geometries.
   void SetTCUseSelectAndMerge(bool use) { fTCUseSelectAndMerge = use; }
   void SetTCUseCircleMerge(bool use) { fTCUseCircleMerge = use; }
   void SetTCCircleMergeTol(double radiusTolFrac, double centerTolMM)
   {
      fTCCircleMergeRadiusTol = radiusTolFrac;
      fTCCircleMergeCenterTol = centerTolMM;
   }

   /// Charge sign from the angular sweep about the fitted circle centre (robust on
   /// shallow arcs). Forwarded to AtPRA::SetChargeFromCenter. Default off.
   void SetChargeFromCenter(bool use) { fChargeFromCenter = use; }

   /// Riemann (PRAlgorithm == 3): RANSAC inlier distance to plane in (x,y,x^2+y^2) space (mm)
   void SetRiemannInlierDist(double d) { fRiemannInlierDist = d; }
   void SetRiemannMinHits(int n) { fRiemannMinHits = n; }
   void SetRiemannMaxTracks(int n) { fRiemannMaxTracks = n; }
   void SetRiemannMaxIter(int n) { fRiemannMaxIter = n; }
   /// Apply the AT-TPC primary/fragment heuristic after Riemann finding.
   /// Disable for annular geometries (PUMA) where the beam-axis vertex assumption is invalid.
   void SetRiemannUseSelectAndMerge(bool use) { fRiemannUseSelectAndMerge = use; }
   void SetRiemannUseArcWalkExtend(bool use) { fRiemannUseArcWalkExtend = use; }
   void SetRiemannArcWalkWindow(int n) { fRiemannArcWalkWindow = n; }
   void SetRiemannArcWalkMaxMiss(int n) { fRiemannArcWalkMaxMiss = n; }

   ClassDef(AtPRAtask, 1);
};

#endif
