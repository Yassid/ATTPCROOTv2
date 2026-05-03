/// @file AtPiGunGenerator.h
/// @brief Random π+/π− gun for UKF mass-discrimination studies.
///
/// Per event: draws a random multiplicity uniformly in [fMinMult, fMaxMult],
/// then emits that many pions, each with random sign (π+ or π−), isotropic
/// direction, and momentum magnitude uniform in [fMinP_GeV, fMaxP_GeV].
/// Vertex geometry mirrors AtPUMAGenerator (Gaussian (x,y) clipped to
/// fTrapRadius cylinder, uniform z in [fMeanZ ± fDeltaZ]) so the same
/// digi+reco chain can run unmodified.
///
/// Intended use: isolate the K-vs-π mass discrimination in the UKF by
/// generating only pions and only fitting pion hypotheses — if species ID
/// is now perfect, mass confusion was driving the PUMA results; if not,
/// the energy-loss / dE/dx model needs work.

#ifndef ATPIGUNGENERATOR_H
#define ATPIGUNGENERATOR_H

#include <FairGenerator.h>
#include <Rtypes.h>

class FairPrimaryGenerator;

class AtPiGunGenerator : public FairGenerator {
public:
   AtPiGunGenerator() = default;
   ~AtPiGunGenerator() override = default;

   Bool_t ReadEvent(FairPrimaryGenerator *primGen) override;

   /// Multiplicity range (inclusive). Default 1..3.
   void SetMultiplicityRange(Int_t lo, Int_t hi)
   {
      fMinMult = lo;
      fMaxMult = hi;
   }
   /// Momentum magnitude range in GeV/c (uniform). Default 0.10..0.30 — covers
   /// the PUMA pion KE spectrum reasonably.
   void SetMomentumRange(Double_t pLo_GeV, Double_t pHi_GeV)
   {
      fMinP_GeV = pLo_GeV;
      fMaxP_GeV = pHi_GeV;
   }
   /// Force a single charge for ALL pions in the event when set: +1 → only
   /// π+, −1 → only π−, 0 (default) → random per pion.
   void SetCharge(Int_t q) { fForcedCharge = q; }

   /// Vertex distribution (mm). Same conventions as AtPUMAGenerator.
   void SetVertexXY(Double_t meanX, Double_t sigmaX, Double_t meanY, Double_t sigmaY)
   {
      fMeanX = meanX;
      fSigmaX = sigmaX;
      fMeanY = meanY;
      fSigmaY = sigmaY;
   }
   void SetVertexZHalfRange(Double_t deltaZ_mm) { fDeltaZ = deltaZ_mm; }
   void SetVertexZOffset(Double_t meanZ_mm) { fMeanZ = meanZ_mm; }
   void SetTrapRadius(Double_t r_mm) { fTrapRadius = r_mm; }

private:
   Int_t fMinMult{1};
   Int_t fMaxMult{3};
   Double_t fMinP_GeV{0.10};
   Double_t fMaxP_GeV{0.30};
   Int_t fForcedCharge{0};

   Double_t fMeanX{0.}, fSigmaX{0.}, fMeanY{0.}, fSigmaY{0.};
   Double_t fDeltaZ{0.};
   Double_t fMeanZ{0.};
   Double_t fTrapRadius{10.};

   ClassDefOverride(AtPiGunGenerator, 1);
};

#endif
