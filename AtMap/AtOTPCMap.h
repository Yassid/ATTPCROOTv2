#ifndef ATOTPCMAP_H
#define ATOTPCMAP_H

#include "AtMap.h"

#include <Math/Point2Dfwd.h>
#include <Rtypes.h>

#include <unordered_map>

class TBuffer;
class TClass;
class TMemberInspector;

class AtOTPCMap : public AtMap {

public:
   AtOTPCMap();
   ~AtOTPCMap();

   void Dump() override;                                         // pure virtual member
   void GeneratePadPlane() override;                             // pure virtual member
   ROOT::Math::XYPoint CalcPadCenter(Int_t PadRef) override;     // pure virtual member
   Int_t BinToPad(Int_t binval) override { return binval - 1; }; // pure virtual member

private:
   std::unordered_map<Int_t, Int_t> fBinToPadTable;
   std::unordered_map<Int_t, Int_t>::iterator fBinToPadTableIt;

   ClassDefOverride(AtOTPCMap, 1);
};

#endif