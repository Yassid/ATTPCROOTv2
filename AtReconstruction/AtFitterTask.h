/*********************************************************************
 *   Fitter Task AtFitterTask.hh			             *
 *   Author: Y. Ayyad ayyadlim@frib.msu.edu            	             *
 *   Log: 3/10/2021 					             *
 *								     *
 *********************************************************************/

#ifndef ATFITTERTASK
#define ATFITTERTASK



#include <FairTask.h>

#include <Rtypes.h>

#include <cstddef>
#include <TClonesArray.h>  // for TClonesArray
#include <TString.h>       // for TString
#include <memory>          // for unique_ptr
#include "AtFitter.h"      // for AtFitter
class AtDigiPar;
class TBuffer;
class TClass;
class TMemberInspector;


class AtFitterTask : public FairTask {

public:
   // AtFitterTask();
   ~AtFitterTask() = default;
   AtFitterTask(std::unique_ptr<AtFITTER::AtFitter> fitter);

   void SetInputBranch(TString branchName);
   void SetOutputBranch(TString branchName);
   void SetPersistence(Bool_t value = kTRUE);

   virtual InitStatus Init();
   virtual void SetParContainers();
   virtual void Exec(Option_t *opt);


private:
   TString fInputBranchName;
   TString fOutputBranchName;

   Bool_t fIsPersistence; //!< Persistence check variable

   std::unique_ptr<AtFITTER::AtFitter> fFitter;
   AtDigiPar *fPar{nullptr};
   TClonesArray *fPatternEventArray;
   TClonesArray fTrackingEventArray;

   std::size_t fEventCnt{0};

   ClassDef(AtFitterTask, 1);
};

#endif
