/*********************************************************************
 *  AtMergerHDFUnpacker  -- see header for the format description.
 *********************************************************************/

#include "AtMergerHDFUnpacker.h"

#include "AtRawEvent.h"

#include <FairLogger.h>

#include <Rtypes.h>
#include <TString.h>

#include <H5Fpublic.h> // H5Fopen / H5Fclose
#include <H5Gpublic.h> // H5Giterate
#include <H5Lpublic.h> // H5Lexists
#include <H5Ppublic.h> // H5P_DEFAULT
#include <H5public.h>  // htri_t

#include <algorithm> // for max_element, min_element
#include <regex>
#include <string>
#include <tuple>
#include <vector>

ClassImp(AtMergerHDFUnpacker);

AtMergerHDFUnpacker::AtMergerHDFUnpacker(mapPtr map) : AtHDFUnpacker(std::move(map)) {}

bool AtMergerHDFUnpacker::IsMergerFile(char const *file)
{
   hid_t f = H5Fopen(file, H5F_ACC_RDONLY, H5P_DEFAULT);
   if (f < 0)
      return false;
   htri_t hasEvents = H5Lexists(f, "events", H5P_DEFAULT);
   htri_t hasGet = H5Lexists(f, "get", H5P_DEFAULT);
   H5Fclose(f);
   // raw merger layout: /events present, legacy /get absent
   return (hasEvents > 0) && (hasGet <= 0);
}

std::size_t AtMergerHDFUnpacker::open(char const *file)
{
   auto f = open_file(file, AtHDFUnpacker::IO_MODE::READ);
   if (f == 0)
      return 0;
   _file = f;

   auto group_n_entries = open_group(f, "events"); // raw layout: pad data lives under /events
   if (std::get<0>(group_n_entries) == -1)
      return 0;
   _group = std::get<0>(group_n_entries);
   setFirstAndLastEventNum();
   return std::get<1>(group_n_entries);
}

void AtMergerHDFUnpacker::setFirstAndLastEventNum()
{
   // The raw layout has no `meta` dataset; event ids are the group names
   // event_<N> under /events. Collect them and take min/max.
   auto collect = [](hid_t, const char *name, void *op_data) -> herr_t {
      std::string digits = std::regex_replace(std::string(name), std::regex("[^0-9]"), "");
      if (!digits.empty())
         static_cast<std::vector<long> *>(op_data)->push_back(std::stol(digits));
      return 0;
   };
   int idx = 0;
   std::vector<long> ids;
   H5Giterate(_file, "events", &idx, collect, &ids);

   if (ids.empty()) {
      LOG(error) << "AtMergerHDFUnpacker: no event_<N> groups found under /events";
      fFirstEvent = 0;
      fLastEvent = 0;
      return;
   }
   fFirstEvent = *std::min_element(ids.begin(), ids.end());
   fLastEvent = *std::max_element(ids.begin(), ids.end());
   LOG(info) << "AtMergerHDFUnpacker: events " << fFirstEvent << " to " << fLastEvent << " (" << ids.size()
             << " present)";
}

void AtMergerHDFUnpacker::setEventIDAndTimestamps()
{
   // The raw GET stream carries no evt<N>_header dataset (timestamps live in the
   // FRIB coincidence block, which the reco does not need). Set none.
   fRawEvent->SetEventID(fEventID);
   fRawEvent->SetNumberOfTimestamps(0);
}

void AtMergerHDFUnpacker::processData()
{
   // Pad data at /events/event_<N>/get_traces, opened relative to _group (=/events).
   TString dataset = TString::Format("event_%lld/get_traces", static_cast<long long>(fDataEventID));
   std::size_t npads = n_pads(dataset.Data());

   for (std::size_t ipad = 0; ipad < npads; ++ipad)
      processPad(ipad);

   if (npads > 0)   // only a successfully-opened dataset needs closing (skip missing events)
      end_raw_event();
}
