#include "AtPIDStampTask.h"

#include "AtFittedTrack.h"
#include "AtPatternEvent.h"
#include "AtTrack.h"
#include "AtTrackingEvent.h"

#include <FairLogger.h>
#include <FairRootManager.h>

#include <TClonesArray.h>

#include <cmath>

AtPIDStampTask::AtPIDStampTask() : FairTask("AtPIDStampTask") {}

void AtPIDStampTask::AddFittedBranch(const TString &n)
{
   if (!fBranchesUserSet) { // first user call replaces the defaults
      fFittedBranches.clear();
      fBranchesUserSet = true;
   }
   fFittedBranches.push_back(n);
}

void AtPIDStampTask::UsePionKaon()
{
   fPID.AddSpecies("pi", 139.57039, 211);
   fPID.AddSpecies("K", 493.677, 321);
}

InitStatus AtPIDStampTask::Init()
{
   FairRootManager *ioMan = FairRootManager::Instance();
   if (ioMan == nullptr) {
      LOG(error) << "AtPIDStampTask: FairRootManager not instantiated!";
      return kERROR;
   }

   fPatternArray = dynamic_cast<TClonesArray *>(ioMan->GetObject(fPatternBranch));
   if (fPatternArray == nullptr) {
      LOG(error) << "AtPIDStampTask: cannot find pattern branch " << fPatternBranch;
      return kERROR;
   }

   fFittedArrays.clear();
   for (const auto &bn : fFittedBranches) {
      auto *arr = dynamic_cast<TClonesArray *>(ioMan->GetObject(bn));
      if (arr == nullptr) {
         LOG(warn) << "AtPIDStampTask: fitted branch " << bn << " not found, skipping";
         continue;
      }
      fFittedArrays.push_back(arr);
      LOG(info) << "AtPIDStampTask: will stamp PID onto " << bn;
   }
   if (fFittedArrays.empty()) {
      LOG(error) << "AtPIDStampTask: no fitted-track branches found to stamp";
      return kERROR;
   }
   if (fPID.GetSpecies().empty()) {
      LOG(warn) << "AtPIDStampTask: no species registered, defaulting to pi + K";
      UsePionKaon();
   }
   LOG(info) << "AtPIDStampTask initialized (B=" << fBField << " T, k=" << fPID.GetCalibration() << ", "
             << fPID.GetSpecies().size() << " species)";
   return kSUCCESS;
}

void AtPIDStampTask::Exec(Option_t * /*opt*/)
{
   if (fPatternArray == nullptr || fPatternArray->GetEntriesFast() == 0)
      return;
   auto *pat = dynamic_cast<AtPatternEvent *>(fPatternArray->At(0));
   if (pat == nullptr)
      return;
   auto &tracks = pat->GetTrackCand(); // std::vector<AtTrack>

   for (auto *arr : fFittedArrays) {
      if (arr == nullptr || arr->GetEntriesFast() == 0)
         continue;
      auto *te = dynamic_cast<AtTrackingEvent *>(arr->At(0));
      if (te == nullptr)
         continue;

      for (const auto &ft : te->GetFittedTracks()) { // const vector, but ft-> is non-const
         int tid = ft->GetTrackID();
         if (tid < 0 || tid >= static_cast<int>(tracks.size()))
            continue;

         double rig = 0, dedx = 0;
         if (!AtTools::AtBetheBlochPID::Observables(tracks[tid], fBField, rig, dedx))
            continue;

         // rigidity -> full momentum for betagamma: p = p_T / sin(theta)
         if (fPitchCorrect) {
            double th = ft->GetKinematics(0).theta;
            double s = std::sin(th);
            if (th > 0 && std::abs(s) > 0.2)
               rig = rig / std::abs(s);
         }

         int idx = fPID.Classify(rig, dedx);
         if (idx < 0)
            continue;
         const auto &sp = fPID.GetSpecies()[idx];

         // charge sign: prefer the PRA curvature sign (most reliable), else the fit's.
         int sign = tracks[tid].GetChargeSign();
         if (sign == 0)
            sign = (ft->GetParticleInfo(0).charge >= 0) ? 1 : -1;
         std::string name = sp.name + (sign > 0 ? "+" : "-");

         ft->SetParticleInfo(0, name, sign, sp.massMeV);
         ++fStamped;
      }
   }
}

void AtPIDStampTask::Finish()
{
   LOG(info) << "AtPIDStampTask: stamped PID on " << fStamped << " fitted tracks";
}

ClassImp(AtPIDStampTask);
