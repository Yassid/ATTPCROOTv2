#ifndef ATFITTRACKMETADATA_H
#define ATFITTRACKMETADATA_H

#include <Rtypes.h> // for Double_t, THashConsistencyHolder, ClassDefOverride
#include <TObject.h>

class TBuffer;
class TClass;
class TMemberInspector;

/**
 * Class for storing the result of the fit of an AtTrack from an AtFitter class.
 */
class AtFitTrackMetadata : public TObject {
protected:
   // Statistics parameters of the fit.
   Double_t fPValue{0};
   Double_t fChi2{0};
   Int_t fNdf{0};
   Bool_t fFitConverged{false};

   // Post-fit quality flag: kTRUE if the fit passed the quality cut (converged,
   // sensible chi2/ndf, physical kinematics). Tracks that fail are KEPT (not dropped)
   // so curlers/blobs stay inspectable; downstream code can filter on this.
   Bool_t fGoodFit{false};

   // Material-effects provenance of THIS result. A fitter that retries a failed
   // material-effects fit with material effects OFF produces a track from a different
   // model, and mixing the two populations in one spectrum degrades the resolution
   // while looking like a single sample. These two flags make the split explicit:
   //   fMatEffects         - material effects were active in the fit that produced this
   //                         result (kFALSE both when they were never requested and when
   //                         the requested fit fell back).
   //   fMatEffectsFallback - material effects were REQUESTED, that fit threw, and this
   //                         result comes from the no-material retry. Only ever kTRUE in
   //                         a material-effects production, where it isolates exactly the
   //                         tracks to exclude before quoting a resolution.
   Bool_t fMatEffects{false};
   Bool_t fMatEffectsFallback{false};

   // The track ID for which this fit was done for.
   Int_t fTrackID{-1};

public:
   AtFitTrackMetadata() = default;
   ~AtFitTrackMetadata() = default;

   void SetPValue(Double_t value) { fPValue = value; }
   void SetChi2(Double_t value) { fChi2 = value; }
   void SetNdf(Int_t value) { fNdf = value; }
   void SetFitConverged(Bool_t value) { fFitConverged = value; }
   void SetGoodFit(Bool_t value) { fGoodFit = value; }
   void SetMatEffects(Bool_t value) { fMatEffects = value; }
   void SetMatEffectsFallback(Bool_t value) { fMatEffectsFallback = value; }
   void SetTrackID(Int_t value) { fTrackID = value; }

   Double_t GetPValue() const { return fPValue; }
   Double_t GetChi2() const { return fChi2; }
   Int_t GetNdf() const { return fNdf; }
   Bool_t GetFitConverged() const { return fFitConverged; }
   Bool_t GetGoodFit() const { return fGoodFit; }
   Bool_t GetMatEffects() const { return fMatEffects; }
   Bool_t GetMatEffectsFallback() const { return fMatEffectsFallback; }
   Int_t GetTrackID() const { return fTrackID; }

   virtual void Print() const;

   ClassDefOverride(AtFitTrackMetadata, 3);
};

#endif
