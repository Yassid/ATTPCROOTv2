#ifndef AtPSAMAX_H
#define AtPSAMAX_H

#include "AtPSA.h"

#include <Rtypes.h> // for Bool_t, THashConsistencyHolder, ClassDefOverride

#include <array>  // for array
#include <memory> // for make_unique, unique_ptr

class AtPad;
class TBuffer;
class TClass;
class TMemberInspector;

/**
 * @brief Simple max finding PSA method.
 *
 *
 *
 */
class AtPSAMax : public AtPSA {

private:
   Bool_t fIsTimeCorr{false};
   // Subtract this many TBs from the peak index before computing z. Accounts
   // for the AGET response peak lagging the electron impulse by ~u_peak·τ_peak
   // (u_peak≈1.167 for AtNominalResponse). Without it, every reconstructed
   // pad-hit z is biased toward the pad plane by τ_peak/TBTime·u_peak · TBlen.
   // Default 0 (no correction) to preserve legacy behaviour.
   Double_t fPeakingShiftTBs{0.0};

public:
   virtual HitVector AnalyzePad(AtPad *pad) override;
   std::unique_ptr<AtPSA> Clone() override { return std::make_unique<AtPSAMax>(*this); }

   void SetTimeCorrection(Bool_t value) { fIsTimeCorr = value; }
   void SetPeakingShiftTBs(Double_t tbs) { fPeakingShiftTBs = tbs; }

private:
   bool shouldSaveHit(double charge, double threshold, int tb);
   Double_t getTBCorr(std::array<Double_t, 512> &trace, int maxAdcIdx);

   ClassDefOverride(AtPSAMax, 1)
};

#endif
