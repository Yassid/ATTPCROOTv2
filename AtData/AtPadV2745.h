#ifndef _ATPADV2745_H
#define _ATPADV2745_H

#include "AtPadBase.h"

#include <Rtypes.h>

#include <cstddef>
#include <functional> // IWYU pragma: keep
#include <iosfwd>
#include <memory>

class TBuffer;
class TClass;
class TMemberInspector;

struct AtV2745Block {
   UShort_t channel;
   ULong64_t time_ps;
   ULong64_t coarse_time_int;
   UShort_t fine_time_int;
   UShort_t board;
   UShort_t charge;
   UShort_t energy_short; // Check!
   UInt_t flags;
   UShort_t pads;
   Double_t charge_cal;
   UShort_t row;
   UShort_t column;
   UShort_t section;
};
// bool operator<(const AtV2745Block &l, const AtV2745Block &r);//TODO
bool operator==(const AtV2745Block &l, const AtV2745Block &r);
std::ostream &operator<<(std::ostream &os, const AtV2745Block &t);

class AtPadV2745 : public AtPadBase {
private:
   AtV2745Block fRef;

public:
   AtPadV2745(AtV2745Block ref = {}) : fRef(ref) {}
   virtual std::unique_ptr<AtPadBase> Clone() const override { return std::make_unique<AtPadV2745>(*this); }
   AtV2745Block GetReference() { return fRef; }
   ClassDefOverride(AtPadV2745, 1);
};
#endif //#ifndef _PADV2745_H
