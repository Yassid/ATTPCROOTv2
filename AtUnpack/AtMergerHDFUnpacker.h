/*********************************************************************
 *  AtMergerHDFUnpacker
 *
 *  Reads the RAW libattpc_merger "events" HDF5 layout, i.e. files whose
 *  top-level group is `/events` with per-event subgroups:
 *
 *      /events/event_<N>/get_traces          Dataset {npads, 517}  (native short)
 *      /events/event_<N>/frib_physics/1903   FRIB IC traces  (not used here)
 *      /events/event_<N>/frib_physics/977    FRIB coincidence (not used here)
 *      /scalers                              (not used here)
 *
 *  This is the format produced directly by the FRIB DAQ / attpc_merger before
 *  the "remerge" step that flattens everything into the legacy
 *  `/get/evt<N>_data` + `/frib` + `/meta` layout read by AtHDFUnpacker.
 *  (Spyral calls this the MergerV1 format; see spyral trace_reader.py.)
 *
 *  The GET pad row (5 electronics-id words + 512 ADC samples = 517 shorts) is
 *  BYTE-IDENTICAL to the legacy `evt<N>_data`, so all pad decoding is inherited
 *  unchanged from AtHDFUnpacker; only the dataset navigation is overridden.
 *
 *  Backward compatibility: AtHDFUnpacker is not modified. Pick the reader with
 *  AtMergerHDFUnpacker::IsMergerFile(file) (true for the raw layout).
 *********************************************************************/

#ifndef _AtMERGERHDFUNPACKER_H_
#define _AtMERGERHDFUNPACKER_H_

#include "AtHDFUnpacker.h"
#include "AtUnpacker.h" // for mapPtr

#include <Rtypes.h> // for ClassDefOverride

#include <cstddef> // for size_t

class TBuffer;
class TClass;
class TMemberInspector;

class AtMergerHDFUnpacker : public AtHDFUnpacker {

public:
   AtMergerHDFUnpacker(mapPtr map);
   ~AtMergerHDFUnpacker() = default;

   /// @brief True if `file` is a raw libattpc_merger "events" layout:
   ///        it has a top-level `/events` group and no legacy `/get` group.
   static bool IsMergerFile(char const *file);

protected:
   std::size_t open(char const *file) override;
   void setFirstAndLastEventNum() override;
   void processData() override;
   void setEventIDAndTimestamps() override;

   ClassDefOverride(AtMergerHDFUnpacker, 1);
};

#endif
