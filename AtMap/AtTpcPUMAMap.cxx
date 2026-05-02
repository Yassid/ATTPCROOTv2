#include "AtTpcPUMAMap.h"

#include <FairLogger.h>
#include <Math/Point2D.h>
#include <TH2Poly.h>
#include <TMath.h>

#include <cmath>
#include <iostream>

ClassImp(AtTpcPUMAMap);

AtTpcPUMAMap::AtTpcPUMAMap(double rInner_mm, double rOuter_mm, int nRings, int nPadsPerRing)
   : fRIn(rInner_mm), fROut(rOuter_mm), fNRings(nRings), fNPadsPerRing(nPadsPerRing)
{
   const int n = fNRings * fNPadsPerRing;
   AtPadCoord.resize(boost::extents[n][4][2]); // 4 vertices per trapezoidal pad
   std::fill(AtPadCoord.data(), AtPadCoord.data() + AtPadCoord.num_elements(), 0.);
   fNumberPads = n;

   // Equal-area ring radii: r_{i+1}^2 = r_i^2 + (rOut^2 - rIn^2) / nRings
   fRingRadii.resize(fNRings + 1);
   const double r2_step = (fROut * fROut - fRIn * fRIn) / fNRings;
   fRingRadii[0] = fRIn;
   for (int i = 1; i <= fNRings; ++i)
      fRingRadii[i] = std::sqrt(fRingRadii[i - 1] * fRingRadii[i - 1] + r2_step);

   // Populate vertex coordinates so AtPadCoord users (event display etc.) work.
   const double dPhi = 2. * TMath::Pi() / fNPadsPerRing;
   for (int ring = 0; ring < fNRings; ++ring) {
      const double rIn_i = fRingRadii[ring];
      const double rOut_i = fRingRadii[ring + 1];
      for (int sec = 0; sec < fNPadsPerRing; ++sec) {
         const int pad = ring * fNPadsPerRing + sec;
         const double phi0 = (sec - 0.5) * dPhi;
         const double phi1 = (sec + 0.5) * dPhi;
         AtPadCoord[pad][0][0] = rIn_i * std::cos(phi0);
         AtPadCoord[pad][0][1] = rIn_i * std::sin(phi0);
         AtPadCoord[pad][1][0] = rIn_i * std::cos(phi1);
         AtPadCoord[pad][1][1] = rIn_i * std::sin(phi1);
         AtPadCoord[pad][2][0] = rOut_i * std::cos(phi1);
         AtPadCoord[pad][2][1] = rOut_i * std::sin(phi1);
         AtPadCoord[pad][3][0] = rOut_i * std::cos(phi0);
         AtPadCoord[pad][3][1] = rOut_i * std::sin(phi0);
      }
   }
   kIsParsed = true;
}

void AtTpcPUMAMap::Dump()
{
   std::cout << "AtTpcPUMAMap: " << fNRings << " rings x " << fNPadsPerRing << " pads = " << fNumberPads
             << " pads, R in [" << fRIn << ", " << fROut << "] mm." << std::endl;
}

void AtTpcPUMAMap::GeneratePadPlane()
{
   if (fPadPlane) {
      LOG(info) << "AtTpcPUMAMap: pad plane already generated, skipping.";
      return;
   }

   fPadPlane = new TH2Poly();
   fPadPlane->SetName("ATTPC_PUMAPadPlane");
   fPadPlane->SetTitle("PUMA annular pad plane");
   fPadPlane->SetFloat();

   for (int pad = 0; pad < (int)fNumberPads; ++pad) {
      Double_t px[5] = {AtPadCoord[pad][0][0], AtPadCoord[pad][1][0], AtPadCoord[pad][2][0], AtPadCoord[pad][3][0],
                        AtPadCoord[pad][0][0]};
      Double_t py[5] = {AtPadCoord[pad][0][1], AtPadCoord[pad][1][1], AtPadCoord[pad][2][1], AtPadCoord[pad][3][1],
                        AtPadCoord[pad][0][1]};
      fPadPlane->AddBin(5, px, py);
   }
   fPadPlane->ChangePartition(std::max(50, fNPadsPerRing / 4), std::max(50, fNPadsPerRing / 4));
   LOG(info) << "AtTpcPUMAMap: generated " << fNumberPads << " pads (" << fNRings << " rings x " << fNPadsPerRing
             << "), R=[" << fRIn << "," << fROut << "] mm.";
}

ROOT::Math::XYPoint AtTpcPUMAMap::CalcPadCenter(Int_t pad)
{
   if (pad < 0 || pad >= (int)fNumberPads)
      return {-9999., -9999.};
   const int ring = pad / fNPadsPerRing;
   const int sec = pad % fNPadsPerRing;
   const double rMid = 0.5 * (fRingRadii[ring] + fRingRadii[ring + 1]);
   const double phi = sec * (2. * TMath::Pi() / fNPadsPerRing);
   return {rMid * std::cos(phi), rMid * std::sin(phi)};
}
