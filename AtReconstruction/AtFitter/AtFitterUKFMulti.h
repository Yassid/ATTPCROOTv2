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
   /// and is stored in the resulting AtFittedTrack particle info.
   void AddHypothesis(std::unique_ptr<AtFitterUKF> fitter, std::string name)
   {
      fHypotheses.push_back({std::move(fitter), std::move(name)});
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
   };
   std::vector<Hypo> fHypotheses; //! transient, not persisted by ROOT
   double fChi2Cap{1e30};

   ClassDefOverride(AtFitterUKFMulti, 1);
};

} // namespace EventFit

#endif
