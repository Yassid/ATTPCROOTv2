#include "AtGainMatchTask.h"

#include "AtPIDEvent.h"

#include <FairLogger.h>
#include <FairRootManager.h>

#include <TString.h>

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

ClassImp(AtGainMatchTask);

Int_t AtGainMatchTask::RunNumberFromName(const TString &name)
{
   // Take the LAST digit group: "run_0017_reco" -> 17, but also "a1975/run_0017" -> 17.
   // A trailing "_reco"/"_pid" has no digits, so it does not interfere; a name that is all
   // text yields -1 and the caller is expected to treat that as an error.
   const std::string s = name.Data();
   int best = -1;
   for (size_t i = 0; i < s.size();) {
      if (std::isdigit(static_cast<unsigned char>(s[i]))) {
         size_t j = i;
         while (j < s.size() && std::isdigit(static_cast<unsigned char>(s[j])))
            ++j;
         best = std::atoi(s.substr(i, j - i).c_str());
         i = j;
      } else {
         ++i;
      }
   }
   return best;
}

InitStatus AtGainMatchTask::Init()
{
   if (fRun < 0) {
      LOG(error) << "AtGainMatchTask: no run number set. Call SetRun() (or "
                    "SetRun(AtGainMatchTask::RunNumberFromName(fileName))). Refusing to guess.";
      return kERROR;
   }
   if (fTablePath.Length() == 0) {
      LOG(error) << "AtGainMatchTask: no gain-match table set. Call SetTable().";
      return kERROR;
   }

   std::ifstream in(fTablePath.Data());
   if (!in) {
      LOG(error) << "AtGainMatchTask: cannot open gain-match table " << fTablePath;
      return kERROR;
   }
   std::string line;
   int nrows = 0;
   while (std::getline(in, line)) {
      if (line.empty() || line[0] == '#')
         continue;
      if (line.rfind("run,", 0) == 0) // header
         continue;
      std::stringstream ss(line);
      std::string runStr, facStr, srcStr;
      if (!std::getline(ss, runStr, ','))
         continue;
      if (!std::getline(ss, facStr, ','))
         continue;
      std::getline(ss, srcStr, ',');
      const int run = std::atoi(runStr.c_str());
      const double fac = std::atof(facStr.c_str());
      if (!(fac > 0.0)) {
         LOG(error) << "AtGainMatchTask: non-positive factor " << fac << " for run " << run << " in " << fTablePath
                    << ". A negative or zero gain factor would flip or erase every dE/dx.";
         return kERROR;
      }
      fTable[run] = fac;
      fSource[run] = srcStr.empty() ? std::string("?") : srcStr;
      ++nrows;
   }
   if (nrows == 0) {
      LOG(error) << "AtGainMatchTask: gain-match table " << fTablePath << " has no usable rows.";
      return kERROR;
   }

   auto it = fTable.find(fRun);
   if (it == fTable.end()) {
      if (!fAllowMissing) {
         LOG(error) << "AtGainMatchTask: run " << fRun << " is not in " << fTablePath << " (" << nrows
                    << " rows, runs " << fTable.begin()->first << "-" << fTable.rbegin()->first
                    << "). Regenerate the table over a range that covers it, or call "
                       "SetAllowMissingRun(kTRUE) to run UNMATCHED on purpose.";
         return kERROR;
      }
      fFactor = 1.0;
      LOG(warn) << "AtGainMatchTask: run " << fRun << " absent from the table; continuing with factor 1.0 "
                << "-- this run is NOT gain matched and must not be mixed with runs that are.";
   } else {
      fFactor = it->second;
   }

   auto *ioMan = FairRootManager::Instance();
   if (ioMan == nullptr) {
      LOG(error) << "AtGainMatchTask: no FairRootManager.";
      return kERROR;
   }
   fPIDArray = dynamic_cast<TClonesArray *>(ioMan->GetObject(fPIDBranchName));
   if (fPIDArray == nullptr) {
      LOG(error) << "AtGainMatchTask: branch " << fPIDBranchName
                 << " not found. This task must be added AFTER AtPIDTask.";
      return kERROR;
   }

   LOG(info) << "AtGainMatchTask: run " << fRun << " factor " << fFactor << " (" << fSource[fRun] << ") from "
             << fTablePath << "; dEdx *= f, sqrtdEdx *= sqrt(f)" << (fScaleDE ? ", dE *= f" : ", dE UNCHANGED")
             << (fScaleClassic ? ", classic estimator also scaled" : ", Spyral estimator only");
   return kSUCCESS;
}

void AtGainMatchTask::Exec(Option_t *)
{
   ++fNEvents;
   if (fPIDArray == nullptr || fPIDArray->GetEntriesFast() == 0)
      return;

   const double sq = std::sqrt(fFactor);

   for (int i = 0; i < fPIDArray->GetEntriesFast(); ++i) {
      auto *ev = dynamic_cast<AtPIDEvent *>(fPIDArray->At(i));
      if (ev == nullptr)
         continue;

      // AtPIDEvent exposes only const getters, so scale copies and refill through the public
      // Clear/Add API. This keeps the framework class untouched; the vectors are a handful of
      // tracks per event, so the copy is irrelevant next to the fitting that follows.
      auto spyral = ev->GetSpyral();
      auto classic = ev->GetClassic();

      for (auto &r : spyral) {
         r.dEdx *= fFactor;
         r.sqrtdEdx *= sq;
         if (fScaleDE)
            r.dE *= fFactor;
      }
      if (fScaleClassic) {
         for (auto &r : classic) {
            r.dEdx *= fFactor;
            r.sqrtdEdx *= sq;
            if (fScaleDE) {
               r.dE *= fFactor;
               r.elossMean *= fFactor;
               r.elossTrunc *= fFactor;
            }
         }
      }

      ev->Clear();
      for (const auto &r : classic)
         ev->AddClassic(r);
      for (const auto &r : spyral)
         ev->AddSpyral(r);

      fNScaled += static_cast<Long64_t>(spyral.size());
   }
}

void AtGainMatchTask::Finish()
{
   LOG(info) << "AtGainMatchTask: run " << fRun << ", factor " << fFactor << " -- scaled " << fNScaled
             << " Spyral PID tracks over " << fNEvents << " events.";
}
