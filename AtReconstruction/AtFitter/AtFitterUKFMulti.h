/// @file AtFitterUKFMulti.h
/// @brief Multi-particle-hypothesis wrapper around AtFitterUKF.
///
/// Holds an ordered list of pre-configured AtFitterUKF instances (one per
/// particle hypothesis: K+, pi-, p, ...). For each AtTrack, runs every
/// hypothesis through its underlying UKF and returns the AtFittedTrack with
/// the best (smallest) reduced chi^2 / ndf. The chosen hypothesis is recorded
/// in the AtFittedTrack's ParticleInfo.
///
/// Usage:
///   auto multi = std::make_unique<EventFit::AtFitterUKFMulti>();
///   multi->AddHypothesis(makeUKF_Kplus(),  "K+");
///   multi->AddHypothesis(makeUKF_PiMinus(),"pi-");
///   AtFitterTask *task = new AtFitterTask(std::move(multi));

#ifndef ATFITTERUKFMULTI_H
#define ATFITTERUKFMULTI_H

#include "AtFitter.h"
#include "AtFitterUKF.h"

#include <Rtypes.h>

#include <memory>
#include <string>
#include <vector>

namespace EventFit {

class AtFitterUKFMulti : public AtFitter {
public:
   AtFitterUKFMulti() = default;
   ~AtFitterUKFMulti() override = default;

   /// Append a fully-configured fitter as one hypothesis. Name appears in logs
   /// and is stored in the resulting AtFittedTrack particle info. Optional
   /// chargeSign (+1, -1) lets the multi-fitter skip this hypothesis whenever
   /// the input track has a definite charge sign (AtTrack::GetChargeSign() set
   /// by the PRA's curvature) and the two disagree. Pass 0 to always evaluate
   /// the hypothesis (default — preserves the legacy 4-hypothesis behaviour).
   void AddHypothesis(std::unique_ptr<AtFitterUKF> fitter, std::string name, int chargeSign = 0)
   {
      fHypotheses.push_back({std::move(fitter), std::move(name), chargeSign});
   }

   void Init() override;

   /// Optional: minimum acceptable reduced chi^2/ndf — fits below this are
   /// considered valid. Hypotheses with chi^2/ndf above the cap are still
   /// considered if no other converged. Default 1e30 (no cap).
   void SetChi2Cap(double cap) { fChi2Cap = cap; }

   /// Number of hypotheses currently registered.
   std::size_t GetNHypotheses() const { return fHypotheses.size(); }

protected:
   AtFittedTrack *GetFittedTrack(AtTrack *track, AtFitMetadata *fitMetadata = nullptr, AtRawEvent *rawEvent = nullptr,
                                 AtEvent *event = nullptr) override;

private:
   struct Hypo {
      std::unique_ptr<AtFitterUKF> fitter;
      std::string name;
      int chargeSign{0};
   };
   std::vector<Hypo> fHypotheses; //! transient, not persisted by ROOT
   double fChi2Cap{1e30};

   ClassDefOverride(AtFitterUKFMulti, 1);
};

} // namespace EventFit

#endif
