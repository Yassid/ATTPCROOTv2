#include "AtPIDTask.h"

#include "AtHit.h"
#include "AtPIDEvent.h"
#include "AtPIDGainTable.h"
#include "AtPatternEvent.h"
#include "AtTrack.h"

#include <FairLogger.h>
#include <FairRootManager.h>

AtPIDTask::AtPIDTask()
   : FairTask("AtPIDTask"), fPIDArray("AtPIDEvent", 1), fClassicEst(fBField, fSmallPadRadius)
{
   fSpyralEst.SetBField(fBField);
   fHelix3DEst.SetBField(fBField);
}

void AtPIDTask::SetBField(Double_t b)
{
   fBField = b;
   fClassicEst.SetBField(b);
   fSpyralEst.SetBField(b);
   fHelix3DEst.SetBField(b);
}

void AtPIDTask::SetSmallPadRadius(Double_t r)
{
   fSmallPadRadius = r;
   fClassicEst.SetSmallPadRadius(r);
   fSpyralEst.SetSmallPadRadius(r);
   fHelix3DEst.SetSmallPadRadius(r);
}

Bool_t AtPIDTask::SetPadGainMap(TString file, Int_t column)
{
   fPadGain = AtTools::AtPIDGainTable::LoadPadGains(std::string(file.Data()), column);
   if (fPadGain.empty()) {
      LOG(error) << "AtPIDTask: pad-gain map " << file << " not loaded -- pad correction OFF";
      return kFALSE;
   }
   LOG(info) << "AtPIDTask: pad-gain map " << file << " (column " << column << "), " << fPadGain.size() << " pads";
   return kTRUE;
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
      // pad-gain correction on a COPY: the pattern event stays untouched for the fitters downstream
      AtTrack corrected;
      if (!fPadGain.empty()) {
         corrected = track;
         for (auto &h : corrected.GetHitArray()) {
            auto it = fPadGain.find(h->GetPadNum());
            if (it != fPadGain.end())
               h->SetCharge(h->GetCharge() / it->second);
         }
      }
      AtTrack &tr = fPadGain.empty() ? const_cast<AtTrack &>(track) : corrected;
      auto cl = fClassicEst.Estimate(tr);
      cl.trackID = tr.GetTrackID();
      pidEvent->AddClassic(cl);
      // AtSpyralPID stamps trackID/nClusters itself, before its own early returns, so even
      // tracks with an invalid PID come back identified.
      auto sp = fSpyralEst.Estimate(tr);
      pidEvent->AddSpyral(sp);
      // the Spyral direction is one of AtHelix3D's three votes
      if (fRunHelix3D) {
         auto h3 = fHelix3DEst.Estimate(tr, sp.valid ? sp.direction : -1);
         pidEvent->AddHelix3D(h3);
         if (fCombine)
            pidEvent->AddCombined(AtTools::CombinePID(sp, h3, fCombineOpt));
      }
   }
}

ClassImp(AtPIDTask);
