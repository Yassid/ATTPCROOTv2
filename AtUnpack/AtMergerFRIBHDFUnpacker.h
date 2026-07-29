/*********************************************************************
 *  AtMergerFRIBHDFUnpacker
 *
 *  Reads the FRIB auxiliary traces (ion chamber, coincidence register) out of
 *  the RAW libattpc_merger "events" HDF5 layout, i.e. the same files that
 *  AtMergerHDFUnpacker reads for the pad data.
 *
 *  The FRIB payload is present in BOTH formats and is byte-identical; only the
 *  path differs:
 *
 *      legacy remerged :  /frib/evt/evt<N>_1903            Dataset {2048, 8}
 *      raw merger      :  /events/event_<N>/frib_physics/1903   Dataset {2048, 8}
 *                         /events/event_<N>/frib_physics/977    Dataset {1}
 *
 *  Module 1903 holds the 8 generic auxiliary channels, 2048 samples each; the
 *  ION CHAMBER is channel 0 (trace[0]), used for the beam-species gate. Because
 *  the dataset shape and content are identical, ALL decoding is inherited from
 *  AtFRIBHDFUnpacker unchanged -- only the dataset navigation is overridden.
 *
 *  Without this class the 11 raw-merger runs of a2091 have no ion chamber and so
 *  cannot be beam-gated, which would leave the PID plane mixing every beam
 *  species. AtFRIBHDFUnpacker is not modified (zero backcompat risk); pick the
 *  reader with AtMergerHDFUnpacker::IsMergerFile(file).
 *********************************************************************/

#ifndef _AtMERGERFRIBHDFUNPACKER_H_
#define _AtMERGERFRIBHDFUNPACKER_H_

#include "AtFRIBHDFUnpacker.h"
#include "AtUnpacker.h" // for mapPtr

#include <Rtypes.h> // for ClassDefOverride

#include <cstddef> // for size_t

class TBuffer;
class TClass;
class TMemberInspector;

class AtMergerFRIBHDFUnpacker : public AtFRIBHDFUnpacker {

public:
   AtMergerFRIBHDFUnpacker(mapPtr map);
   ~AtMergerFRIBHDFUnpacker() = default;

   void Init() override;

protected:
   std::size_t open(char const *file) override;
   void setFirstAndLastEventNum() override;
   void processData() override;

   ClassDefOverride(AtMergerFRIBHDFUnpacker, 1);
};

#endif
