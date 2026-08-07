/// @file AtBeamHole.h
/// @brief Disable the pads inside the AT-TPC beam hole in a SIMULATION.
///
/// WHY THIS IS NEEDED. The AT-TPC pad plane has a beam hole; the region is not read out. A
/// simulation, though, digitizes the beam in full and fills pads the real detector never sees.
/// The consequence is not cosmetic: the beam point cloud sits right on top of the reaction
/// vertex, so the pattern recognition merges the beam into the recoil-particle cluster. On
/// a1954 14C(p,p') that was 67 % of all recoil protons swallowed by the beam cluster, which in
/// turn destroyed the PID (a merged cluster fits a near-straight line, so Brho blows up) and
/// the fit. Removing the hole pads fixed it: MERGED 67.1 % -> 0.2 %, CLEAN 28.8 % -> 97.4 %,
/// Spyral <radius> 250 mm -> 72 mm against 85 mm in the data.
///
/// HOW. Pads whose centre lies inside holeR are inhibited as AtMap::InhibitType::kTotal, which
/// AtPulse honours (AtPulse.cxx:92) -- so they never produce a signal at all, rather than being
/// filtered downstream. Call it AFTER GeneratePadPlane() and BEFORE building AtPulseTask.
///
///   #include "../AtBeamHole.h"
///   auto mapping = std::make_shared<AtTpcMap>();
///   mapping->ParseXMLMap(mapFile);
///   mapping->GeneratePadPlane();
///   AtSim::InhibitBeamHole(mapping, 30.0);   // a1954: 3 cm
///
/// CHOOSING holeR -- it is per experiment, do not copy 30 mm blindly. Measure it from the data:
/// histogram the hit radius for a data run and for the sim, normalise each to its own total,
/// and take the radius beyond which the sim/data ratio goes flat. On a1954 run_0056 the ratio
/// was 4260x at r < 4 mm, 40x at 16-20 mm, and flat (~0.36) past ~30 mm -> holeR = 30 mm,
/// which matched the known 3 cm hole. **15C and 16C (a1975) still need this and have not been
/// checked** -- their sims currently digitize the beam over the hole just as 14C did.
///
/// Returns the number of pads inhibited (0 means something is wrong -- see the warning).

#ifndef ATBEAMHOLE_H
#define ATBEAMHOLE_H

// Interpreted macros get AtMap from the ROOT dictionary autoload; an explicit include would
// require ROOT_INCLUDE_PATH to contain build/include, which build/config.sh does NOT set. The
// sim macros run fine without it today, so including unconditionally here would break them.
#ifndef __CLING__
#include "AtMap.h"
#endif

#include <cmath>
#include <iostream>
#include <memory>

namespace AtSim {

inline int InhibitBeamHole(std::shared_ptr<AtMap> map, double holeR, bool verbose = true)
{
   if (!map) {
      std::cout << "\033[1;31mInhibitBeamHole: null map\033[0m" << std::endl;
      return 0;
   }
   if (holeR <= 0) {
      if (verbose)
         std::cout << "BEAM HOLE: disabled (holeR = " << holeR << ")" << std::endl;
      return 0;
   }

   const unsigned nPads = map->GetNumPads();
   if (nPads == 0) {
      std::cout << "\033[1;31mInhibitBeamHole: map has 0 pads -- call GeneratePadPlane() first"
                << "\033[0m" << std::endl;
      return 0;
   }

   int nOff = 0;
   for (unsigned padNum = 0; padNum < nPads; ++padNum) {
      auto c = map->CalcPadCenter(padNum); // mm
      if (std::isnan(c.X()) || std::isnan(c.Y()))
         continue;
      if (std::hypot(c.X(), c.Y()) < holeR) {
         map->InhibitPad(padNum, AtMap::InhibitType::kTotal);
         ++nOff;
      }
   }

   if (verbose)
      std::cout << "BEAM HOLE: inhibited " << nOff << " of " << nPads << " pads with r < " << holeR << " mm"
                << std::endl;
   // Loud, not silent: a hole that quietly inhibits nothing looks exactly like a working one
   // until the clustering results come out wrong days later.
   if (nOff == 0)
      std::cout << "\033[1;31mWARNING: beam hole requested (r < " << holeR
                << " mm) but NO pad was inhibited -- check CalcPadCenter/the map\033[0m" << std::endl;

   return nOff;
}

} // namespace AtSim

#endif
