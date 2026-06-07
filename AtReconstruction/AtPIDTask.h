#ifndef ATPIDTASK_H
#define ATPIDTASK_H

#include "AtPIDEstimator.h"
#include "AtSpyralPID.h"

#include <FairTask.h>

#include <Rtypes.h>
#include <TClonesArray.h>
#include <TString.h>

/**
 * @brief Runs the two PID estimators on every pattern-recognized track.
 *
 * Reads the AtPatternEvent branch, and for each AtTrack computes PID observables
 * with BOTH the legacy AtTools::AtPIDEstimator (charge/arclength) and the faithful
 * Spyral-style AtTools::AtSpyralPID (first-arc + spline arclength). Results are
 * written to an AtPIDEvent branch (parallel per-track vectors) so the two methods
 * can be compared with full statistics. The legacy estimator is kept untouched.
 */
class AtPIDTask : public FairTask {
private:
   TClonesArray *fPatternEventArray{nullptr}; ///< input (AtPatternEvent)
   TClonesArray fPIDArray;                     ///< output (AtPIDEvent)

   TString fInputBranchName{"AtPatternEvent"};
   TString fOutputBranchName{"AtPIDEvent"};
   Bool_t fIsPersistence{kTRUE};

   Double_t fBField{2.85};        ///< Tesla
   Double_t fSmallPadRadius{152.0};

   AtTools::AtPIDEstimator fClassicEst;
   AtTools::AtSpyralPID fSpyralEst;

public:
   AtPIDTask();
   virtual ~AtPIDTask() = default;

   void SetBField(Double_t b);
   void SetSmallPadRadius(Double_t r);
   void SetPersistence(Bool_t v) { fIsPersistence = v; }
   void SetInputBranch(TString n) { fInputBranchName = n; }
   void SetOutputBranch(TString n) { fOutputBranchName = n; }

   InitStatus Init() override;
   void Exec(Option_t *opt) override;

   ClassDefOverride(AtPIDTask, 1);
};

#endif // ATPIDTASK_H
