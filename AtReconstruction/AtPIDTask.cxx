#include "AtPIDTask.h"

#include "AtPIDEvent.h"
#include "AtPatternEvent.h"
#include "AtTrack.h"

#include <FairLogger.h>
#include <FairRootManager.h>

AtPIDTask::AtPIDTask()
   : FairTask("AtPIDTask"), fPIDArray("AtPIDEvent", 1), fClassicEst(fBField, fSmallPadRadius)
{
   fSpyralEst.SetBField(fBField);
}

void AtPIDTask::SetBField(Double_t b)
{
   fBField = b;
   fClassicEst.SetBField(b);
   fSpyralEst.SetBField(b);
}

void AtPIDTask::SetSmallPadRadius(Double_t r)
{
   fSmallPadRadius = r;
   fClassicEst.SetSmallPadRadius(r);
   fSpyralEst.SetSmallPadRadius(r);
}

InitStatus AtPIDTask::Init()
{
   FairRootManager *ioMan = FairRootManager::Instance();
   if (ioMan == nullptr) {
      LOG(error) << "AtPIDTask: FairRootManager not instantiated!";
      return kERROR;
   }

   fPatternEventArray = dynamic_cast<TClonesArray *>(ioMan->GetObject(fInputBranchName));
   if (fPatternEventArray == nullptr) {
      LOG(error) << "AtPIDTask: cannot find input branch " << fInputBranchName;
      return kERROR;
   }

   ioMan->Register(fOutputBranchName, "AtTPC", &fPIDArray, fIsPersistence);
   LOG(info) << "AtPIDTask initialized: " << fInputBranchName << " -> " << fOutputBranchName << " (B=" << fBField
             << " T)";
   return kSUCCESS;
}

void AtPIDTask::Exec(Option_t * /*opt*/)
{
   fPIDArray.Delete();
   if (fPatternEventArray->GetEntriesFast() == 0)
      return;

   auto *patternEvent = dynamic_cast<AtPatternEvent *>(fPatternEventArray->At(0));
   if (patternEvent == nullptr)
      return;

   auto *pidEvent = new (fPIDArray[0]) AtPIDEvent();
   // One entry per pattern track, in order, in BOTH vectors -- downstream code may rely on that
   // parallelism. It should still prefer matching on trackID, because a fitted track records only
   // its ID and AtTrackFinderTC's ID is a cluster label, not an array position.
   for (auto &track : patternEvent->GetTrackCand()) {
      AtTrack &tr = const_cast<AtTrack &>(track);
      auto cl = fClassicEst.Estimate(tr);
      cl.trackID = tr.GetTrackID();
      pidEvent->AddClassic(cl);
      // AtSpyralPID stamps trackID/nClusters itself, before its own early returns, so even
      // tracks with an invalid PID come back identified.
      pidEvent->AddSpyral(fSpyralEst.Estimate(tr));
   }
}

ClassImp(AtPIDTask);
