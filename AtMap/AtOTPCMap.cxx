#include "AtOTPCMap.h"

#include <FairLogger.h>

#include <Math/Point2D.h>
#include <Rtypes.h>
#include <TH2Poly.h>
#include <TString.h>

#include <boost/multi_array/base.hpp>
#include <boost/multi_array/extent_gen.hpp>
#include <boost/multi_array/multi_array_ref.hpp>
#include <boost/multi_array/subarray.hpp>

#include <algorithm>
#include <cstdlib>
#include <fstream> // IWYU pragma: keep
#include <iostream>
#include <string>
#include <vector>

using XYPoint = ROOT::Math::XYPoint;

AtOTPCMap::AtOTPCMap() : AtMap()
{
   AtPadCoord.resize(boost::extents[1024][4][2]);
   std::fill(AtPadCoord.data(), AtPadCoord.data() + AtPadCoord.num_elements(), 0);
   std::cout << " OTPC Map initialized " << std::endl;
   std::cout << " OTPC Pad Coordinates container initialized " << std::endl;
   fNumberPads = 1012;
}

AtOTPCMap::~AtOTPCMap() = default;

void AtOTPCMap::Dump() {}

void AtOTPCMap::GeneratePadPlane()
{
   if (fPadPlane) {
      LOG(error) << "Skipping generation of pad plane, it is already parsed!";
      return;
   }

   Float_t pad_size = 2.0;      // mm
   Float_t pad_spacing = 0.001; // mm

   std::vector<int> pads_per_row{12, 12, 12, 11, 11, 10, 10, 9, 8, 7, 6, 3};

   int pad_num = 0;

   for (auto irow = 0; irow < pads_per_row.size(); ++irow) {

      for (auto ipad = 0; ipad < pads_per_row[irow]; ++ipad) {
         AtPadCoord[ipad + pad_num][0][0] = ipad * pad_size;
         AtPadCoord[ipad + pad_num][0][1] = -irow * pad_size;
         AtPadCoord[ipad + pad_num][1][0] = ipad * pad_size + pad_size;
         AtPadCoord[ipad + pad_num][1][1] = -irow * pad_size;
         AtPadCoord[ipad + pad_num][2][0] = ipad * pad_size + pad_size;
         AtPadCoord[ipad + pad_num][2][1] = -pad_size - irow * pad_size;
         AtPadCoord[ipad + pad_num][3][0] = ipad * pad_size;
         AtPadCoord[ipad + pad_num][3][1] = -pad_size - irow * pad_size;
      }

      std::cout << " Row " << irow << " Number of pads " << pad_num << "\n";
      pad_num += pads_per_row[irow];
   }

   for (auto irow = 0; irow < pads_per_row.size(); ++irow) {

      for (auto ipad = 0; ipad < pads_per_row[irow]; ++ipad) {
         AtPadCoord[ipad + pad_num][0][0] = -ipad * pad_size;
         AtPadCoord[ipad + pad_num][0][1] = -irow * pad_size;
         AtPadCoord[ipad + pad_num][1][0] = -ipad * pad_size - pad_size;
         AtPadCoord[ipad + pad_num][1][1] = -irow * pad_size;
         AtPadCoord[ipad + pad_num][2][0] = -ipad * pad_size - pad_size;
         AtPadCoord[ipad + pad_num][2][1] = -pad_size - irow * pad_size;
         AtPadCoord[ipad + pad_num][3][0] = -ipad * pad_size;
         AtPadCoord[ipad + pad_num][3][1] = -pad_size - irow * pad_size;
      }

      std::cout << " Row " << irow << " Number of pads " << pad_num << "\n";
      pad_num += pads_per_row[irow];
   }

   for (auto irow = 0; irow < pads_per_row.size(); ++irow) {

      for (auto ipad = 0; ipad < pads_per_row[irow]; ++ipad) {
         AtPadCoord[ipad + pad_num][0][0] = -ipad * pad_size;
         AtPadCoord[ipad + pad_num][0][1] = irow * pad_size;
         AtPadCoord[ipad + pad_num][1][0] = -ipad * pad_size - pad_size;
         AtPadCoord[ipad + pad_num][1][1] = irow * pad_size;
         AtPadCoord[ipad + pad_num][2][0] = -ipad * pad_size - pad_size;
         AtPadCoord[ipad + pad_num][2][1] = pad_size + irow * pad_size;
         AtPadCoord[ipad + pad_num][3][0] = -ipad * pad_size;
         AtPadCoord[ipad + pad_num][3][1] = pad_size + irow * pad_size;
      }

      std::cout << " Row " << irow << " Number of pads " << pad_num << "\n";
      pad_num += pads_per_row[irow];
   }

   pad_num += 1;

   for (auto irow = 0; irow < pads_per_row.size(); ++irow) {

      for (auto ipad = 0; ipad < pads_per_row[irow]; ++ipad) {
         AtPadCoord[ipad + pad_num][0][0] = ipad * pad_size;
         AtPadCoord[ipad + pad_num][0][1] = irow * pad_size;
         AtPadCoord[ipad + pad_num][1][0] = ipad * pad_size + pad_size;
         AtPadCoord[ipad + pad_num][1][1] = irow * pad_size;
         AtPadCoord[ipad + pad_num][2][0] = ipad * pad_size + pad_size;
         AtPadCoord[ipad + pad_num][2][1] = pad_size + irow * pad_size;
         AtPadCoord[ipad + pad_num][3][0] = ipad * pad_size;
         AtPadCoord[ipad + pad_num][3][1] = pad_size + irow * pad_size;
      }

      std::cout << " Row " << irow << " Number of pads " << pad_num << "\n";
      pad_num += pads_per_row[irow];
   }

   //////////////////////////////////////////////////////////////////////////////////////////

   std::cout << " Total pads " << pad_num << "\n";

   // fPadInd = pad_num;
   kIsParsed = true;

   fPadPlane = new TH2Poly(); // NOLINT
   for (auto ipad = 0; ipad < pad_num; ++ipad) {
      Double_t px[] = {AtPadCoord[ipad][0][0], AtPadCoord[ipad][1][0], AtPadCoord[ipad][2][0], AtPadCoord[ipad][3][0],
                       AtPadCoord[ipad][0][0]};
      Double_t py[] = {AtPadCoord[ipad][0][1], AtPadCoord[ipad][1][1], AtPadCoord[ipad][2][1], AtPadCoord[ipad][3][1],
                       AtPadCoord[ipad][0][1]};
      fPadPlane->AddBin(5, px, py);
   }

   fPadPlane->SetName("OTPC_Plane");
   fPadPlane->SetTitle("OTPC_Plane");
   fPadPlane->ChangePartition(500, 500);
}

XYPoint AtOTPCMap::CalcPadCenter(Int_t PadRef)
{

   if (!kIsParsed) {
      LOG(error) << " AtOTPCMap::CalcPadCenter Error : Pad plane has not been generated or parsed";
      return {-9999, -9999};
   }

   if (PadRef == -1) { // Boost multi_array crashes with a negative index
      LOG(debug) << " AtOTPCMap::CalcPadCenter Error : Pad not found";
      return {-9999, -9999};
   }

   Float_t x = (AtPadCoord[PadRef][0][0] + AtPadCoord[PadRef][1][0]) / 2.0;
   Float_t y = (AtPadCoord[PadRef][1][1] + AtPadCoord[PadRef][2][1]) / 2.0;
   return {x, y};
}