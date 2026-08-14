#include "AtTpcSquareMap.h"

#include <FairLogger.h>
#include <Math/Point2D.h>
#include <TH2Poly.h>

#include <cmath>
#include <iostream>
#include <limits>

ClassImp(AtTpcSquareMap);

AtTpcSquareMap::AtTpcSquareMap(double padSize_mm, int nPadsX, int nPadsY,
                               double originX_mm, double originY_mm)
   : fPadSize_mm(padSize_mm), fNPadsX(nPadsX), fNPadsY(nPadsY),
     fHalfX(0.5 * nPadsX * padSize_mm), fHalfY(0.5 * nPadsY * padSize_mm),
     fOriginX(std::isnan(originX_mm) ? -fHalfX : originX_mm),
     fOriginY(std::isnan(originY_mm) ? -fHalfY : originY_mm)
{
   const int n = fNPadsX * fNPadsY;
   AtPadCoord.resize(boost::extents[n][4][2]); // 4 vertices per square
   std::fill(AtPadCoord.data(), AtPadCoord.data() + AtPadCoord.num_elements(), 0.);
   fNumberPads = n;

   // Populate vertex coordinates so AtPadCoord users (e.g. event display) work.
   for (int iy = 0; iy < fNPadsY; ++iy) {
      for (int ix = 0; ix < fNPadsX; ++ix) {
         const int pad = ix + iy * fNPadsX;
         const double x0 = fOriginX + ix * fPadSize_mm;
         const double y0 = fOriginY + iy * fPadSize_mm;
         AtPadCoord[pad][0][0] = x0;
         AtPadCoord[pad][0][1] = y0;
         AtPadCoord[pad][1][0] = x0 + fPadSize_mm;
         AtPadCoord[pad][1][1] = y0;
         AtPadCoord[pad][2][0] = x0 + fPadSize_mm;
         AtPadCoord[pad][2][1] = y0 + fPadSize_mm;
         AtPadCoord[pad][3][0] = x0;
         AtPadCoord[pad][3][1] = y0 + fPadSize_mm;
      }
   }
   // Synthetic electronics references, one per pad.
   //
   // WITHOUT THIS THE MAP IS UNUSABLE WITH ANY PAD-LEVEL BOOKKEEPING. AtMap::IsInhibited and
   // InhibitPad both go through GetPadRef(padNum), which reads fPadMapInverse; a map that never
   // fills it returns the SAME default-constructed AtPadReference for every pad, so all 62500
   // pads share one key. Inhibiting a single pad -- as the beam-hole helper does -- then marks
   // the whole plane as kTotal, AtPulse::GetGain returns 0 for every electron, and the
   // digitisation silently produces an empty event ("Skipped 100% of N points"). Measured before
   // this was added: 100 % of points skipped with a 20 mm beam hole, 0 % with the hole disabled.
   //
   // The values are a bookkeeping key only -- this map has no real electronics behind it -- but
   // they are laid out in the AT-TPC's own 68 channels / 4 agets / 4 asads nesting so that
   // anything reasoning about the tuple sees plausible ranges.
   for (int pad = 0; pad < n; ++pad) {
      AtPadReference ref;
      ref.ch = pad % 68;
      ref.aget = (pad / 68) % 4;
      ref.asad = (pad / (68 * 4)) % 4;
      ref.cobo = pad / (68 * 4 * 4);
      fPadMap[ref] = pad;
      fPadMapInverse[pad] = ref;
   }

   kIsParsed = true;
}

void AtTpcSquareMap::Dump()
{
   std::cout << "AtTpcSquareMap: " << fNPadsX << " x " << fNPadsY << " = " << fNumberPads << " pads, "
             << fPadSize_mm << " mm pitch, active " << 2 * fHalfX << " x " << 2 * fHalfY << " mm^2." << std::endl;
}

void AtTpcSquareMap::GeneratePadPlane()
{
   if (fPadPlane) {
      LOG(info) << "AtTpcSquareMap: pad plane already generated, skipping.";
      return;
   }

   fPadPlane = new TH2Poly();
   fPadPlane->SetName("ATTPC_SquarePadPlane");
   fPadPlane->SetTitle("Square pad plane");
   fPadPlane->SetFloat();

   for (int pad = 0; pad < (int)fNumberPads; ++pad) {
      // Close the polygon: 5 points = 4 vertices + return to first.
      Double_t px[5] = {AtPadCoord[pad][0][0], AtPadCoord[pad][1][0], AtPadCoord[pad][2][0], AtPadCoord[pad][3][0],
                        AtPadCoord[pad][0][0]};
      Double_t py[5] = {AtPadCoord[pad][0][1], AtPadCoord[pad][1][1], AtPadCoord[pad][2][1], AtPadCoord[pad][3][1],
                        AtPadCoord[pad][0][1]};
      fPadPlane->AddBin(5, px, py);
   }
   fPadPlane->ChangePartition(std::max(50, fNPadsX), std::max(50, fNPadsY));
   LOG(info) << "AtTpcSquareMap: generated " << fNumberPads << " square pads (" << fPadSize_mm << " mm pitch).";
}

ROOT::Math::XYPoint AtTpcSquareMap::CalcPadCenter(Int_t pad)
{
   if (pad < 0 || pad >= (int)fNumberPads)
      return {-9999., -9999.};
   const int ix = pad % fNPadsX;
   const int iy = pad / fNPadsX;
   const double cx = fOriginX + (ix + 0.5) * fPadSize_mm;
   const double cy = fOriginY + (iy + 0.5) * fPadSize_mm;
   return {cx, cy};
}
