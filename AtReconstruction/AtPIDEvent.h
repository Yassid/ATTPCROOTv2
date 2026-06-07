#ifndef ATPIDEVENT_H
#define ATPIDEVENT_H

#include "AtPIDEstimator.h" // for AtTools::AtPIDResult
#include "AtSpyralPID.h"    // for AtTools::AtSpyralResult

#include <TNamed.h>

#include <vector>

/**
 * @brief Per-event PID container holding observables from BOTH estimators.
 *
 * For every pattern-recognized track, AtPIDTask fills one entry in fClassic (the
 * legacy AtTools::AtPIDEstimator: total charge / arclength) and one in fSpyral
 * (the faithful Spyral-style AtTools::AtSpyralPID: first-arc + spline arclength).
 * The two vectors are parallel (index i = the i-th track), so the two PID methods
 * can be compared track-by-track with full statistics straight from the tree.
 */
class AtPIDEvent : public TNamed {
private:
   std::vector<AtTools::AtPIDResult> fClassic;   ///< legacy AtPIDEstimator, one per track
   std::vector<AtTools::AtSpyralResult> fSpyral; ///< Spyral-style AtSpyralPID, one per track

public:
   AtPIDEvent() : TNamed("AtPIDEvent", "AtPIDEvent") {}

   void Clear(Option_t * = "") override
   {
      fClassic.clear();
      fSpyral.clear();
   }

   void AddClassic(const AtTools::AtPIDResult &r) { fClassic.push_back(r); }
   void AddSpyral(const AtTools::AtSpyralResult &r) { fSpyral.push_back(r); }

   const std::vector<AtTools::AtPIDResult> &GetClassic() const { return fClassic; }
   const std::vector<AtTools::AtSpyralResult> &GetSpyral() const { return fSpyral; }

   /// The adopted DEFAULT PID method is the Spyral-style first-arc + spline-arclength
   /// estimator (AtSpyralPID): it removes the dEdx-inflation artifact of the legacy
   /// charge/arclength method and yields a single clean proton band. The legacy
   /// AtPIDEstimator is kept available via GetClassic() for comparison / fallback.
   const std::vector<AtTools::AtSpyralResult> &GetDefault() const { return fSpyral; }

   ClassDefOverride(AtPIDEvent, 1);
};

#endif // ATPIDEVENT_H
