#include "AtPadV2745.h"

#include <iostream>

std::ostream &operator<<(std::ostream &os, const AtV2745Block &t)
{
   // os << "[" << t.cobo << "," << t.asad << "," << t.aget << "," << t.ch << "]";
   return os;
}

bool operator==(const AtV2745Block &l, const AtV2745Block &r)
{
   return l.channel == r.channel && l.time_ps == r.time_ps && l.coarse_time_int == r.coarse_time_int &&
          l.fine_time_int == r.fine_time_int && l.board == r.board && l.charge == r.charge && l.flags == r.flags &&
          l.pads == r.pads && l.charge_cal == r.charge_cal && l.row == r.row && l.column == r.column &&
          l.section == r.section;
}

ClassImp(AtPadV2745);