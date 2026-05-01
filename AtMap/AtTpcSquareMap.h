/*********************************************************************
 *  AtTpcSquareMap — uniform square-pad map for HELIOS-style sandboxes
 *  No XML / lookup needed: pad layout generated programmatically.
 *********************************************************************/
#ifndef ATTPCSQUAREMAP_H
#define ATTPCSQUAREMAP_H

#include "AtMap.h"

#include <Math/Point2Dfwd.h>
#include <Rtypes.h>

class AtTpcSquareMap : public AtMap {
public:
   /// @param padSize_mm  Side length of each square pad.
   /// @param nPadsX      Number of pads along x.
   /// @param nPadsY      Number of pads along y.
   /// Active area is (nPadsX * padSize) x (nPadsY * padSize), centered at (0,0).
   AtTpcSquareMap(double padSize_mm = 2.0, int nPadsX = 100, int nPadsY = 100);
   ~AtTpcSquareMap() override = default;

   void Dump() override;
   void GeneratePadPlane() override;
   ROOT::Math::XYPoint CalcPadCenter(Int_t PadRef) override;
   Int_t BinToPad(Int_t binval) override { return binval - 1; }

   double GetPadSize() const { return fPadSize_mm; }
   int GetNPadsX() const { return fNPadsX; }
   int GetNPadsY() const { return fNPadsY; }

private:
   double fPadSize_mm;
   int fNPadsX;
   int fNPadsY;
   double fHalfX; // half active extent in x (mm)
   double fHalfY; // half active extent in y (mm)

   ClassDefOverride(AtTpcSquareMap, 1);
};

#endif
