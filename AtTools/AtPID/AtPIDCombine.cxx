#include "AtPIDCombine.h"

#include <cmath>

using namespace AtTools;

AtPIDCombinedResult AtTools::CombinePID(const AtSpyralResult &sp, const AtHelix3DResult &h3, const AtPIDCombineOpt &opt)
{
   AtPIDCombinedResult c;
   c.trackID = sp.trackID >= 0 ? sp.trackID : h3.trackID;
   const double b3 = (opt.useArc && h3.brhoArc > 0) ? h3.brhoArc : h3.brho;
   const bool has3D = h3.valid && b3 > 0 && h3.sqrtdEdx > 0;
   auto use3D = [&]() {
      c.source = 1;
      c.brho = b3;
      c.sqrtdEdx = h3.sqrtdEdx * opt.addBackScale;
   };
   if (!sp.valid || sp.brho <= 0) {
      if (opt.fallback3D && has3D)
         use3D();
      return c;
   }
   c.vetoed = has3D && std::log(sp.brho / b3) >= opt.lnRatioMax && (sp.brho < opt.maxBrho || sp.sqrtdEdx > opt.highSq);
   if (!c.vetoed) {
      c.source = 0;
      c.brho = sp.brho;
      c.sqrtdEdx = sp.sqrtdEdx;
   } else {
      use3D();
   }
   return c;
}
