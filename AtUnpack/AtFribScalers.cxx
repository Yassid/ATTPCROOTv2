#include "AtFribScalers.h"

#include <FairLogger.h>

#include <H5Cpp.h> // for H5File, Group, DataSet

#include <cstdint> // for uint32_t
#include <cstdio>  // for printf

ClassImp(AtFribScalers);

namespace {
// The two layouts a1975-era files come in. Legacy remerged files put the scalers under the frib
// group; raw libattpc_merger files put them at the top level.
const char *const kGroupCandidates[] = {"frib/scaler", "scalers"};
} // namespace

bool AtFribScalers::ReadRun(const std::string &fileName)
{
   // The scaler count is recorded nowhere, so the end is found by scanning until a dataset is
   // missing. That miss is the NORMAL termination, not an error, and without this HDF5 prints a
   // twenty-line diagnostic for every run.
   H5::Exception::dontPrint();

   H5::H5File file;
   try {
      file = H5::H5File(fileName.c_str(), H5F_ACC_RDONLY);
   } catch (...) {
      LOG(warn) << "AtFribScalers: cannot open " << fileName;
      return false;
   }

   H5::Group group;
   bool haveGroup = false;
   for (const auto *cand : kGroupCandidates) {
      try {
         group = file.openGroup(cand);
         haveGroup = true;
         break;
      } catch (...) {
      }
   }
   if (!haveGroup) {
      LOG(warn) << "AtFribScalers: no scaler group in " << fileName;
      return false;
   }

   long nEvt = 0;
   long long sum[kNScalers] = {};
   while (true) {
      uint32_t buf[32] = {0};
      try {
         H5::DataSet ds = group.openDataSet(("scaler" + std::to_string(nEvt) + "_data").c_str());
         ds.read(buf, H5::PredType::NATIVE_UINT32);
      } catch (...) {
         break; // the expected end of the scan
      }
      for (int i = 0; i < kNScalers; ++i)
         sum[i] += static_cast<long long>(buf[i]);
      ++nEvt;
   }
   if (nEvt == 0) {
      LOG(warn) << "AtFribScalers: scaler group present but empty in " << fileName;
      return false;
   }

   for (int i = 0; i < kNScalers; ++i)
      fCounts[i] += sum[i];
   fScalerEvents += nEvt;
   ++fRuns;
   return true;
}

void AtFribScalers::Add(const AtFribScalers &other)
{
   for (int i = 0; i < kNScalers; ++i)
      fCounts[i] += other.fCounts[i];
   fScalerEvents += other.fScalerEvents;
   fRuns += other.fRuns;
}

double AtFribScalers::GetLiveFractionClock() const
{
   return fCounts[kClockFree] > 0 ? static_cast<double>(fCounts[kClockLive]) / fCounts[kClockFree] : 1.0;
}

double AtFribScalers::GetLiveFractionTrigger() const
{
   return fCounts[kTriggerFree] > 0 ? static_cast<double>(fCounts[kTriggerLive]) / fCounts[kTriggerFree] : 1.0;
}

const char *AtFribScalers::GetScalerName(EScaler which)
{
   static const char *const names[kNScalers] = {"clock_free",   "clock_live", "trigger_free", "trigger_live",
                                                "ic_sca",       "mesh_sca",   "si1_cfd",      "si2_cfd",
                                                "sipm",         "ic_ds",      "ic_cfd"};
   return (which >= 0 && which < kNScalers) ? names[which] : "unknown";
}

void AtFribScalers::Print() const
{
   printf("  AtFribScalers: %d runs, %ld scaler events\n", fRuns, fScalerEvents);
   for (int i = 0; i < kNScalers; ++i)
      printf("    %-14s %20lld\n", GetScalerName(static_cast<EScaler>(i)), fCounts[i]);
   printf("    live fraction  clock %.4f   trigger %.4f\n", GetLiveFractionClock(), GetLiveFractionTrigger());
}
