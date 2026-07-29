/*********************************************************************
 *  AtMergerFRIBHDFUnpacker -- see header for the format description.
 *********************************************************************/

#include "AtMergerFRIBHDFUnpacker.h"

#include "AtRawEvent.h"

#include <FairLogger.h>

#include <Rtypes.h>
#include <TString.h>

#include <H5Gpublic.h> // H5Giterate
#include <H5Ipublic.h> // hid_t

#include <algorithm> // for max_element, min_element
#include <regex>
#include <string>
#include <tuple>
#include <vector>

ClassImp(AtMergerFRIBHDFUnpacker);

AtMergerFRIBHDFUnpacker::AtMergerFRIBHDFUnpacker(mapPtr map) : AtFRIBHDFUnpacker(std::move(map)) {}

void AtMergerFRIBHDFUnpacker::Init()
{
   auto numEvents = open(fInputFileName.c_str());
   auto uniqueEvents = GetNumEvents();
   LOG(info) << "AtMergerFRIBHDFUnpacker: " << numEvents << " entries under /events, " << uniqueEvents
             << " unique events";
   if (fEventID > uniqueEvents)
      LOG(fatal) << "Exceeded valid range of event numbers. Looking for " << fEventID << ", max event number is "
                 << uniqueEvents;

   // NOTE: deliberately no "numEvents/3" consistency check here. That check belongs to the
   // legacy layout, where /frib/evt holds three flat datasets per event (evt<N>_1903,
   // evt<N>_977, plus a header) so the entry count is 3x the event count. In the raw merger
   // layout /events holds exactly ONE subgroup per event, so entries == events.

   fDataEventID = fFirstEvent + fEventID;
}

std::size_t AtMergerFRIBHDFUnpacker::open(char const *file)
{
   auto f = open_file(file, AtHDFUnpacker::IO_MODE::READ);
   if (f == 0)
      return 0;
   _file = f;

   // Keep _group pointing at /events; the per-event FRIB dataset is then addressed
   // with the relative path event_<N>/frib_physics/1903 in processData().
   auto group_n_entries = open_group(f, "events");
   if (std::get<0>(group_n_entries) == -1) {
      LOG(error) << "AtMergerFRIBHDFUnpacker: no /events group in " << file;
      return 0;
   }
   _group = std::get<0>(group_n_entries);
   setFirstAndLastEventNum();
   return std::get<1>(group_n_entries);
}

void AtMergerFRIBHDFUnpacker::setFirstAndLastEventNum()
{
   // No /meta group in the raw layout: the event ids are the event_<N> subgroup names.
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
      LOG(error) << "AtMergerFRIBHDFUnpacker: no event_<N> groups found under /events";
      fFirstEvent = 0;
      fLastEvent = 0;
      return;
   }
   fFirstEvent = *std::min_element(ids.begin(), ids.end());
   fLastEvent = *std::max_element(ids.begin(), ids.end());
   LOG(info) << "AtMergerFRIBHDFUnpacker: events " << fFirstEvent << " to " << fLastEvent << " (" << ids.size()
             << " present)";
}

void AtMergerFRIBHDFUnpacker::processData()
{
   // Same module (1903), same {2048, 8} shape as the legacy evt<N>_1903 -- only nested one
   // level deeper. n_pads()/processPad()/pad_raw_data() are inherited untouched.
   TString dsName = TString::Format("event_%lld/frib_physics/1903", fDataEventID);
   fRawEvent->SetEventName(dsName.Data());

   std::size_t nch = n_pads(dsName.Data());
   if (nch == 0) {
      LOG(debug) << "AtMergerFRIBHDFUnpacker: no FRIB data for event " << fDataEventID;
      return;
   }

   for (std::size_t ich = 0; ich < nch; ++ich)
      processPad(ich);

   end_raw_event(); // close dataset
}
