#ifndef AtPRA_H
#define AtPRA_H

#include "AtTrack.h" // for AtTrack
#include "AtTrackTransformer.h"

#include <Rtypes.h>  // for Double_t, Float_t, Int_t, THashConsistencyHolder
#include <TObject.h> // for TObject

#include <algorithm> // for max
#include <memory>
#include <type_traits> // for false_type, is_signed, true_type
#include <vector>      // for vector

class AtDigiPar;
class AtEvent;
class AtHit;
class AtPatternEvent;
class TBuffer;
class TClass;
class TMemberInspector;

namespace AtPATTERN {

/**
 * @brief Find patterns in hit clouds.
 *
 * Base class for finding tracks in a hit cloud. Right now, just supports HC.
 *
 *
 */
class AtPRA : public TObject {
protected:
   std::vector<AtTrack> fTrackCand; //< Candidate tracks

   AtDigiPar *fPar; ///< parameter container

   Int_t fMaxHits{5000};
   Int_t fMinHits{0};
   Float_t fMeanDistance{1e9};
   Int_t fKNN{5};             //<! Number of nearest neighbors kNN
   Double_t fStdDevMulkNN{0}; //<! Std dev multiplier for kNN
   Double_t fkNNDist{10};     //<! Distance threshold for outlier rejection in kNN

   Bool_t kSetPrunning{false}; //<<! Enable prunning of tracks

   std::unique_ptr<AtTools::AtTrackTransformer> fTrackTransformer{std::make_unique<AtTools::AtTrackTransformer>()};
   Double_t fClusterRadius{0};      //<! Radius of hit clusters
   Double_t fClusterDistance{0};    //<! Distance between hit clusters
   Double_t fRadiusFitFraction{1.0}; //<! Max fraction of hits for radius fit
   Int_t fMinHitsRadius{3};         //<! Minimum hits for radius fit (robust circle)
   Int_t fMaxHitsRadius{1000};         //<! Maximum hits for radius fit (avoid spiral)
   Bool_t fUseHitsForRadius{false};  //<! Use raw hits instead of clusters for radius

public:
   virtual ~AtPRA() = default;

   // Getters
   virtual std::vector<AtTrack> GetTrackCand() const { return fTrackCand; }

   // Setters
   void SetMaxHits(Int_t maxHits) { fMaxHits = maxHits; }
   void SetMinHits(Int_t minHits) { fMinHits = minHits; }
   void SetMeanDistance(Float_t meanDistance) { fMeanDistance = meanDistance; }
   void SetkNN(Double_t knn) { fKNN = knn; }
   void SetStdDevMulkNN(Double_t stdDevMul) { fStdDevMulkNN = stdDevMul; }
   void SetkNNDist(Double_t dist) { fkNNDist = dist; }
   void SetPrunning() { kSetPrunning = kTRUE; }
   void SetClusterRadius(Double_t clusterRadius) { fClusterRadius = clusterRadius; }
   void SetClusterDistance(Double_t clusterDistance) { fClusterDistance = clusterDistance; }
   /// Max fraction of hits (from vertex end) used for circle radius fit
   void SetRadiusFitFraction(Double_t frac) { fRadiusFitFraction = frac; }
   /// Min/max hits for radius fit (adaptive: uses min(max, frac*total) but at least min)
   void SetMinHitsRadius(Int_t n) { fMinHitsRadius = n; }
   void SetMaxHitsRadius(Int_t n) { fMaxHitsRadius = n; }
   /// Use raw hits instead of clusters for circle fit (default: false)
   void SetUseHitsForRadius(Bool_t use) { fUseHitsForRadius = use; }
   /// Set diffusion parameters for cluster covariance calculation.
   void SetDiffusionParams(double coefT, double coefL, double driftVel, double tbTime, double padResXY = -1)
   {
      fTrackTransformer->SetDiffusionParams(coefT, coefL, driftVel, tbTime, padResXY);
   }

   virtual std::unique_ptr<AtPatternEvent> FindTracks(AtEvent &event) = 0;

   void PruneTrack(AtTrack &track);
   bool kNN(const std::vector<std::unique_ptr<AtHit>> &hits, AtHit &hit, int k);

   /// Reorder clusters along the track arc using nearest-neighbor walk
   /// from the vertex end (highest Z).
   void OrderClustersAlongTrack(AtTrack &track);

protected:
   // Functions that need to be moved to another class. They assume a curved track
   /*
    * Takes track and sets fGeo... parameters using SampleConsensus. I think these are then used
    * as initial guesses for GenFit. Probably, this should be moved into the fitting classes instead of
    * here. At the very least, it needs to be in a subclass that only deals with curved tracks, it
    * would make no sense to apply this function to straight tracks.
    */
   void SetTrackInitialParameters(AtTrack &track);

   template <typename T>
   inline constexpr int GetSign(T num, std::true_type is_signed)
   {
      return (T(0) < num) - (num < T(0));
   }
   template <typename T>
   inline constexpr int GetSign(T num, std::false_type is_signed)
   {
      return (T(0) < num);
   }
   template <typename T>
   inline constexpr int GetSign(T num)
   {
      return GetSign(num, std::is_signed<T>());
   }

   ClassDef(AtPRA, 1)
};

} // namespace AtPATTERN

Double_t fitf(Double_t *x, Double_t *par);

#endif
