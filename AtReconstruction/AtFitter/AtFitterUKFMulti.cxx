#include "AtFitterUKFMulti.h"

#include "AtFitMetadata.h"
#include "AtFittedTrack.h"
#include "AtFitTrackMetadata.h"
#include "AtTrack.h"

#include <FairLogger.h>

#include <limits>
#include <utility>

ClassImp(EventFit::AtFitterUKFMulti);

namespace EventFit {

void AtFitterUKFMulti::Init()
{
   for (auto &h : fHypotheses)
      if (h.fitter)
         h.fitter->Init();
   LOG(info) << "AtFitterUKFMulti: initialised with " << fHypotheses.size() << " hypotheses";
}

AtFittedTrack *AtFitterUKFMulti::GetFittedTrack(AtTrack *track, AtFitMetadata *fitMetadata, AtRawEvent *rawEvent,
                                                AtEvent *event)
{
   if (fHypotheses.empty()) {
      LOG(warn) << "AtFitterUKFMulti::GetFittedTrack: no hypotheses configured";
      return nullptr;
   }

   AtFittedTrack *best = nullptr;
   double bestChi = std::numeric_limits<double>::infinity();
   std::string bestName;

   // Each hypothesis sees a fresh copy of the track so per-fitter mutations
   // (adaptive re-clustering, geo overrides) don't leak between hypotheses.
   for (auto &h : fHypotheses) {
      AtTrack copy(*track);
      AtFittedTrack *cand = h.fitter->GetFittedTrack(&copy, fitMetadata, rawEvent, event);
      if (cand == nullptr)
         continue;

      double chi = std::numeric_limits<double>::infinity();
      const auto &meta = cand->GetTrackMetadata();
      if (meta) {
         double ndf = meta->GetNdf();
         if (ndf > 0)
            chi = meta->GetChi2() / ndf;
      }

      const bool prefer = (best == nullptr) || (chi < bestChi);
      if (prefer) {
         delete best;
         best = cand;
         bestChi = chi;
         bestName = h.name;
      } else {
         delete cand;
      }
   }

   if (best && !bestName.empty()) {
      // Record which hypothesis won. Charge / mass already set by AtFitterUKF;
      // we just override the PDG-name string for downstream analysis.
      const auto &pi = best->GetParticleInfo();
      best->SetParticleInfo(bestName, pi.charge, pi.mass);
   }

   if (best == nullptr)
      LOG(debug) << "AtFitterUKFMulti: all " << fHypotheses.size() << " hypotheses failed for this track";
   return best;
}

} // namespace EventFit
