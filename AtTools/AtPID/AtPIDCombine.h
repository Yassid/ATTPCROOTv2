#ifndef ATPIDCOMBINE_H
#define ATPIDCOMBINE_H

#include "AtHelix3D.h"
#include "AtSpyralPID.h"

#include <Rtypes.h>

namespace AtTools {

/// One combined PID entry per pattern track (AtPIDEvent::GetCombined, filled when AtPIDTask::SetCombine is on).
struct AtPIDCombinedResult {
   Int_t trackID{-1};
   Int_t source{-1};    ///< 0 = Spyral, 1 = AtHelix3D (vetoed by the Spyral check and re-added), -1 = no PID
   Double_t brho{0};    ///< [T m]
   Double_t sqrtdEdx{0}; ///< on the Spyral scale (the 3D value times AtPIDCombineOpt::addBackScale)
   Bool_t vetoed{false}; ///< the Spyral Brho failed the cross-check against the 3D circle
};

/// Parameters of the combination. Defaults reproduce the a1975 D2 recipe (2026-11-01) except addBackScale and
/// highSq, which are detector/gain specific and default to 1 / off.
struct AtPIDCombineOpt {
   /// Spyral Brho rejected when ln(BrhoSpyral / Brho3D) >= lnRatioMax AND (BrhoSpyral < maxBrho OR Spyral
   /// sqrt(dE/dx) > highSq). Slow, near-90-deg tracks are nearly flat in z, Spyral's z-ordered spline is ill-posed
   /// and its radius comes out ~1.6x too large: elastic deuterons then fake a triton band (a1975: 14 % of tracks).
   double lnRatioMax{0.25};
   double maxBrho{0.45};
   double highSq{1e9};      ///< a1975 D2: 30 (Spyral run+pad-gain units) catches the upper part of the fake streak
   double addBackScale{1.0}; ///< 3D sqrt(dE/dx) -> Spyral scale (a1975 D2: 0.965, median Spyral/3D where Brho agree)
   bool useArc{true};        ///< prefer AtHelix3DResult::brhoArc when it was computed
   bool fallback3D{false};   ///< use the 3D result when Spyral has no valid estimate
};

/// The combination rule, kept free of FairRoot so it can be unit tested.
AtPIDCombinedResult CombinePID(const AtSpyralResult &sp, const AtHelix3DResult &h3, const AtPIDCombineOpt &opt);

} // namespace AtTools
#endif
