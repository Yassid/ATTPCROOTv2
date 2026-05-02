/*********************************************************************
 *  AtTpcPUMAMap — annular pad plane for the PUMA TPC sandbox.
 *
 *  Mirrors puma-tpc-simulation/Drift/PadPlane.cpp:
 *    - Default 4096 trapezoidal pads
 *    - 16 concentric equal-area rings, 256 pads per ring
 *    - Inner radius 62.9 mm, outer 121.1 mm
 *  All four constructor args are configurable.
 *
 *  Pad indexing: pad = ring * fNPadsPerRing + sector,
 *                 ring   in [0, fNRings),
 *                 sector in [0, fNPadsPerRing).
 *********************************************************************/
#ifndef ATTPCPUMAMAP_H
#define ATTPCPUMAMAP_H

#include "AtMap.h"

#include <Math/Point2Dfwd.h>
#include <Rtypes.h>

#include <vector>

class AtTpcPUMAMap : public AtMap {
public:
   /// @param rInner_mm  Inner radius of the innermost ring (mm).
   /// @param rOuter_mm  Outer radius of the outermost ring (mm).
   /// @param nRings     Number of concentric equal-area rings.
   /// @param nPadsPerRing  Azimuthal subdivision (pads per ring).
   AtTpcPUMAMap(double rInner_mm = 62.9, double rOuter_mm = 121.1, int nRings = 16, int nPadsPerRing = 256);
   ~AtTpcPUMAMap() override = default;

   void Dump() override;
   void GeneratePadPlane() override;
   ROOT::Math::XYPoint CalcPadCenter(Int_t PadRef) override;
   Int_t BinToPad(Int_t binval) override { return binval - 1; }

   double GetInnerRadius() const { return fRIn; }
   double GetOuterRadius() const { return fROut; }
   int GetNRings() const { return fNRings; }
   int GetNPadsPerRing() const { return fNPadsPerRing; }

private:
   double fRIn;
   double fROut;
   int fNRings;
   int fNPadsPerRing;

   std::vector<double> fRingRadii; // size = fNRings + 1, ring i spans [fRingRadii[i], fRingRadii[i+1]]

   ClassDefOverride(AtTpcPUMAMap, 1);
};

#endif
