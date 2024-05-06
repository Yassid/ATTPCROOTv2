#include "AtSOLARISUnpacker.h"

#include "AtMap.h"
#include "AtPad.h"
#include "AtPadBase.h" // for AtPadBase
#include "AtPadV2745.h"
#include "AtPadValue.h"
#include "AtRawEvent.h"
#include "AtTpcMap.h"

#include <FairLogger.h>

#include <Rtypes.h>

#include "SolReader.h"

#include <algorithm>
#include <array> // for array
#include <future>
#include <iostream>
#include <iterator> // for begin, end
#include <numeric>  // for accumulate
#include <thread>
#include <utility>

AtSOLARISUnpacker::AtSOLARISUnpacker(mapPtr map, Long64_t numBlocks, std::string fileName)
   : AtUnpacker(map), fNumBlocks(numBlocks), fInputFileName(fileName)
{
   fNumBlocks = numBlocks;
}

void AtSOLARISUnpacker::Init()
{
   fReader = std::unique_ptr<SolReader>(new SolReader(fInputFileName));
}

void AtSOLARISUnpacker::FillRawEvent(AtRawEvent &event)
{
   fRawEvent = &event;
   processData();

   fRawEvent->SetIsGood(kTRUE);
   LOG(debug) << " Unpacked " << fRawEvent->GetNumPads() << " pads";
   fEventID++;
   fDataEventID++;
}

void AtSOLARISUnpacker::SetWindowSize(ULong64_t windowSize)
{
   fWindowSize = windowSize;
}

void AtSOLARISUnpacker::processData()
{

   Double_t Coarse_Time0_micros; // Coarse time stamp in micros (to store last coarse timestamps read in files)
   Coarse_Time0_micros = 1e11;   // Initialization value (> 1 day)
   UInt_t HiPFlags0;             // High priority flags (8 bits used) (Intermediate variable to produce single flags)
   UInt_t LoPFlags0;             // Low priority flags (12 bits used) (Intermediate variable to produce single flags)
   ULong64_t Time0_ps;           // Timestamp (Coarse+Fine) in ps

   // array of variables for sorting the fine time)
   UShort_t array_Channel[100];
   ULong64_t array_Time_ps[100];
   ULong64_t array_Coarse_Time[100];
   UShort_t array_Fine_Time_int[100];
   UShort_t array_Board[100];
   UShort_t array_Charge[100];       // Check!
   UShort_t array_Energy_Short[100]; // Check!
   UInt_t array_Flags[100];
   UShort_t array_Pads[100];
   Double_t array_Charge_cal[100];

   ULong64_t nwritten0 = 0;
   ULong64_t nwritten1 = 0;

   int k = 0; // this index is the number of event with the same coarse gain

   fReader->ReadNextBlock(); // read one binary word (the first one)

   array_Channel[k] = (UShort_t)fReader->hit->channel;
   array_Charge[k] = (UShort_t)fReader->hit->charge;
   array_Energy_Short[k] = (UShort_t)fReader->hit->energy_short;
   array_Coarse_Time[k] = (Double_t)fReader->hit->timestamp;
   // printf(" Timestamp_raw = %15.6f\n", (Double_t) reader->hit->timestamp);
   array_Fine_Time_int[k] = (UShort_t)fReader->hit->fine_timestamp;
   array_Time_ps[k] =
      ULong64_t(8 * ns_to_ps * array_Coarse_Time[k]) + ULong64_t((Double_t)array_Fine_Time_int[k] * fineTS_to_ps);
   LoPFlags0 = (UInt_t)fReader->hit->flags_low_priority;
   HiPFlags0 = (UInt_t)fReader->hit->flags_high_priority;
   array_Flags[k] = LoPFlags0 << 16;
   array_Flags[k] += HiPFlags0;
   // array_Board[k] = Board0; // Check!  //TODO

   //#  charge calibration initialisation #
   calib_init(verbosity);
   //#  initialise pad mapping #
   pad_init(verbosity);

   std::cout << "Charge calibration and pad mapping initialised" << std::endl;

   ULong64_t timeinit = array_Time_ps[k];
   std::cout << " time init: " << timeinit << std::endl;

   //##################### Charge calibration ############################
   for (UInt_t cp = 0; cp < 64; cp++) {
      if (array_Board[k] == board_id[0] && array_Channel[k] == cp) {
         array_Charge_cal[k] = array_Charge[k] * cal_fac[0][cp]; // Charge_cal is the Charge calibration
         array_Pads[k] = pad[cp]; // a tracker pad number is associated to the corresponding digitizer channel
      }
      if (array_Board[k] == board_id[1] && array_Channel[k] == cp) {
         array_Charge_cal[k] = array_Charge[k] * cal_fac[1][cp];
         array_Pads[k] = pad[cp];
      }
      if (array_Board[k] == board_id[2] && array_Channel[k] == cp) {
         array_Charge_cal[k] = array_Charge[k] * cal_fac[2][cp];
         array_Pads[k] = pad[cp];
      }
      if (array_Board[k] == board_id[3] && array_Channel[k] == cp) {
         array_Charge_cal[k] = array_Charge[k] * cal_fac[3][cp];
         array_Pads[k] = pad[cp];
      }
      if (array_Board[k] == board_id[4] && array_Channel[k] == cp) {
         array_Charge_cal[k] = array_Charge[k] * cal_fac[4][cp];
         array_Pads[k] = pad[cp];
      }
   }
   //##################### END Charge calibration ############################

   for (int i = 1; i > 0; i++) {
      if (fReader->IsEndOfFile() == true)
         break; // 2024-03-20 added, as in converter_solaris.C

      // if(i % 10000 == 0) std::cout<<"Reading next block : "<< i <<std::endl;

      fReader->ReadNextBlock(); // 2024-03-20 added, as in converter_solaris.C
      k++;

      array_Channel[k] = (UShort_t)fReader->hit->channel;
      array_Charge[k] = (UShort_t)fReader->hit->charge;
      array_Energy_Short[k] = (UShort_t)fReader->hit->energy_short;
      array_Coarse_Time[k] = (Double_t)fReader->hit->timestamp;
      array_Fine_Time_int[k] = (UShort_t)fReader->hit->fine_timestamp;
      array_Time_ps[k] =
         ULong64_t(8 * ns_to_ps * array_Coarse_Time[k]) + ULong64_t(array_Fine_Time_int[k] * fineTS_to_ps);
      LoPFlags0 = (UInt_t)fReader->hit->flags_low_priority;
      HiPFlags0 = (UInt_t)fReader->hit->flags_high_priority;
      array_Flags[k] = LoPFlags0 << 16;
      array_Flags[k] += HiPFlags0;
      // array_Board[k]=Board0;

      std::cout << " Time ps: " << array_Time_ps[k] << " - Timeinit  : " << timeinit
                << " -  Difference : " << array_Time_ps[k] - timeinit << " - Coarse time " << array_Coarse_Time[k]
                << std::endl;
      if ((array_Time_ps[k] - timeinit) > fWindowSize)
         timeinit = array_Time_ps[k];

      //##################### Charge calibration ############################

      for (UInt_t cp = 0; cp < 64; cp++) {
         if (array_Board[k] == board_id[0] && array_Channel[k] == cp) {
            array_Charge_cal[k] = array_Charge[k] * cal_fac[0][cp]; // Charge_cal is the Charge calibration
            array_Pads[k] = pad[cp]; // a tracker pad number is associated to the corresponding digitizer channel
         }
         if (array_Board[k] == board_id[1] && array_Channel[k] == cp) {
            array_Charge_cal[k] = array_Charge[k] * cal_fac[1][cp];
            array_Pads[k] = pad[cp];
         }
         if (array_Board[k] == board_id[2] && array_Channel[k] == cp) {
            array_Charge_cal[k] = array_Charge[k] * cal_fac[2][cp];
            array_Pads[k] = pad[cp];
         }
         if (array_Board[k] == board_id[3] && array_Channel[k] == cp) {
            array_Charge_cal[k] = array_Charge[k] * cal_fac[3][cp];
            array_Pads[k] = pad[cp];
         }
         if (array_Board[k] == board_id[4] && array_Channel[k] == cp) {
            array_Charge_cal[k] = array_Charge[k] * cal_fac[4][cp];
            array_Pads[k] = pad[cp];
         }
      }
      //##################### END Charge calibration ############################

      if (array_Coarse_Time[k] != array_Coarse_Time[k - 1]) {
         if (k == 1) { // the entry has a coarse time different with that of the previous entry.

            AtV2745Block v2745_block;

            v2745_block.coarse_time_int = (ULong64_t)(array_Coarse_Time[k - 1] * us_to_ps); // 2024-03-20 Check!
            // Single_Coarse_Time_int =  (ULong64_t)(array_Coarse_Time[k-1]*ns_to_ps); // 2024-03-20 Check!
            v2745_block.board = array_Board[k - 1];
            v2745_block.channel = array_Channel[k - 1];
            v2745_block.charge = array_Charge[k - 1];
            v2745_block.fine_time_int = array_Fine_Time_int[k - 1];
            v2745_block.time_ps = array_Time_ps[k - 1];
            v2745_block.flags = array_Flags[k - 1];
            v2745_block.charge_cal = array_Charge_cal[k - 1];
            v2745_block.pads = array_Pads[k - 1];
            v2745_block.row = findRow(v2745_block.board, v2745_block.pads);
            v2745_block.section = findSection(v2745_block.board);

            array_Coarse_Time[k - 1] = array_Coarse_Time[k];
            array_Board[k - 1] = array_Board[k];
            array_Channel[k - 1] = array_Channel[k];
            array_Charge[k - 1] = array_Charge[k];
            array_Fine_Time_int[k - 1] = array_Fine_Time_int[k];
            array_Time_ps[k - 1] = array_Time_ps[k];
            array_Flags[k - 1] = array_Flags[k];

            array_Charge_cal[k - 1] = array_Charge_cal[k];
            array_Pads[k - 1] = array_Pads[k];

            /*auto pad = fRawEvent->AddAuxPad(std::to_string(v2745_block.channel)+"_"+std::to_string(i)).first;
            pad->AddAugment("Block", std::make_unique<AtPadV2745>(v2745_block));*/

            nwritten0++;
            k = 0;

         } else if (k > 1) { // if more than one entry has the same coarse time. Order them by the fine time!
                             // sorting of the elements of the array by fine time
            quicksort(array_Fine_Time_int, array_Coarse_Time, array_Board, array_Channel, array_Charge, array_Time_ps,
                      array_Flags, array_Charge_cal, array_Pads, 0, k - 1);

            for (int h = 0; h < k; h++) {
               AtV2745Block v2745_block;
               v2745_block.channel = array_Channel[h];
               v2745_block.charge = array_Charge[h];
               v2745_block.time_ps = array_Time_ps[h];
               v2745_block.coarse_time_int = (ULong64_t)(array_Coarse_Time[h] * us_to_ps);
               v2745_block.fine_time_int = array_Fine_Time_int[h];
               v2745_block.flags = array_Flags[h];
               v2745_block.board = array_Board[h];
               v2745_block.charge_cal = array_Charge_cal[h];
               v2745_block.pads = array_Pads[h];
               v2745_block.row = findRow(v2745_block.board, v2745_block.pads);
               v2745_block.section = findSection(v2745_block.board);

               /*auto pad = fRawEvent->AddAuxPad(std::to_string(v2745_block.channel)+"_"+std::to_string(i)).first;
               pad->AddAugment("Block", std::make_unique<AtPadV2745>(v2745_block));*/

               nwritten0++;
            }
            // make the last entry of the array the first entry of a new group
            array_Coarse_Time[0] = array_Coarse_Time[k];
            array_Board[0] = array_Board[k];
            array_Channel[0] = array_Channel[k];
            array_Charge[0] = array_Charge[k];
            array_Fine_Time_int[0] = array_Fine_Time_int[k];
            array_Time_ps[0] = array_Time_ps[k];
            array_Flags[0] = array_Flags[k];
            array_Charge_cal[0] = array_Charge_cal[k];
            array_Pads[0] = array_Pads[k];

            k = 0;
         }

      } // Different coarse time
   }
}

void AtSOLARISUnpacker::ProcessFile() {}

bool AtSOLARISUnpacker::IsLastEvent()
{
   return 0;
}

void AtSOLARISUnpacker::SetInputFileName(std::string fileName)
{
   fInputFileName = fileName;
}

Long64_t AtSOLARISUnpacker::GetNumEvents()
{

   return 0;
}

void AtSOLARISUnpacker::ProcessInputFile()
{
   LOG(info) << " Processing input file : " << fInputFileName << "\n";
}

Int_t AtSOLARISUnpacker::findRow(int dig, int pad)
{
   int row = -1;
   if (dig == board_id[0]) {
      if (pad == 1000) {
         row = 5;
      } else {
         row = 0;
      }
   } else if (dig == board_id[1]) {
      if (pad == 1000) {
         row = 6;
      } else {
         row = 1;
      }
   } else if (dig == board_id[2]) {
      if (pad == 1000) {
         row = 7;
      } else {
         row = 2;
      }
   } else if (dig == board_id[3]) {
      if (pad == 1000) {
         row = 8;
      } else {
         row = 3;
      }
   } else if (dig == board_id[4]) {
      if (pad == 1000) {
         row = 9;
      } else if (pad == 1001) {
         row = 10;
      } else {
         row = 4;
      }
   }
   return row;
}

Int_t AtSOLARISUnpacker::findColumn(int dig)
{ // 2024-03-20 Check this function!!!
   int Column = -1;
   if (dig == board_id[0]) {
      Column = 0;
   } else if (dig == board_id[1]) {
      Column = 1;
   } else if (dig == board_id[2]) {
      Column = 2;
   } else if (dig == board_id[3]) {
      Column = 3;
   } else if (dig == board_id[4]) {
      Column = 4;
   }
   return Column;
}

Int_t AtSOLARISUnpacker::findSection(int dig)
{
   int section = -1;
   if (dig == board_id[0]) {
      section = 0;
   } else if (dig == board_id[1]) {
      section = 4;
   } else if (dig == board_id[2]) {
      section = 8;
   } else if (dig == board_id[3]) {
      section = 12;
   } else if (dig == board_id[4]) {
      section = 16;
   }
   return section;
}

void AtSOLARISUnpacker::calib_init(int ver)
{
   std::ifstream chargeCalibFile; // charge calibration parameter file
   chargeCalibFile.open(fileCalCharge);
   char buffer[200]; // buffer for text
   int PA_ID;        // preamplifier (PA) number
   //##  read header ##

   for (int i = 0; i < 8; i++) {
      chargeCalibFile.getline(buffer, 200);
      if (ver > 0)
         std::cout << buffer << std::endl;
   }
   //##  read  parameters ##
   std::cout << " read calibration parameters" << std::endl;
   for (int k = 0; k < 5; k++) {            // loop for 5 PAs
      chargeCalibFile.getline(buffer, 200); // read 1st empty line
      if (ver > 0)
         std::cout << buffer << std::endl;
      chargeCalibFile.getline(buffer, 200); // read 2nd line with PA & Dig. info
      if (ver > 0)
         std::cout << buffer << std::endl;
      chargeCalibFile >> PA_ID; // read 3rd line with PA ID
      if (ver > 0)
         std::cout << PA_ID << std::endl;
      chargeCalibFile.getline(buffer, 200); // read 3rd line space
      chargeCalibFile.getline(buffer, 200); // read 4th empty line
      if (ver > 0)
         std::cout << buffer << std::endl;
      for (int i = 0; i < 64; i++) { // loop for the 64 channels of one PA
         chargeCalibFile >> cal_fac[k][i];
         if (ver > 0)
            std::cout << cal_fac[k][i] << std::endl;
      }
      chargeCalibFile.getline(buffer, 200); // read last line space
   }
   chargeCalibFile.close();
}

void AtSOLARISUnpacker::pad_init(int ver)
{
   std::ifstream fp1;
   fp1.open(filePadMapping);

   char buffer[200]; // buffer for text
   int channel;
   //##  read header ##
   for (int i = 0; i < 9; i++) {
      fp1.getline(buffer, 200);
      if (ver > 0)
         std::cout << buffer << std::endl;
   }
   //##  read parameters ##
   for (int i = 0; i < 64; i++) {
      fp1 >> channel >> pad[i];
      if (ver > 0)
         std::cout << channel << " \t " << pad[i] << std::endl;
   }
   fp1.close();
}

void AtSOLARISUnpacker::quicksort(UShort_t arr0[], ULong64_t arr1[], UShort_t arr2[], UShort_t arr3[], UShort_t arr4[],
                                  ULong64_t arr5[], UInt_t arr6[], Double_t arr7[], UShort_t arr8[], int low, int high)
{
   int i = low;
   int j = high;
   UShort_t y0, y2, y3, y4, y8;
   y0 = y2 = y3 = y4 = y8 = 0;
   ULong64_t y5 = 0;
   ULong64_t y1;
   Double_t y7;
   y1 = y7 = 0;
   UInt_t y6 = 0;
   // compare value
   UShort_t z = arr0[(low + high) / 2];

   // partition
   do {
      /* find member above ... */
      while (arr0[i] < z)
         i++;

      /* find element below ... */
      while (arr0[j] > z)
         j--;

      if (i <= j) {
         /* swap two elements */
         y0 = arr0[i];
         arr0[i] = arr0[j];
         arr0[j] = y0;

         y1 = arr1[i];
         arr1[i] = arr1[j];
         arr1[j] = y1;

         y2 = arr2[i];
         arr2[i] = arr2[j];
         arr2[j] = y2;

         y3 = arr3[i];
         arr3[i] = arr3[j];
         arr3[j] = y3;

         y4 = arr4[i];
         arr4[i] = arr4[j];
         arr4[j] = y4;

         y5 = arr5[i];
         arr5[i] = arr5[j];
         arr5[j] = y5;

         y6 = arr6[i];
         arr6[i] = arr6[j];
         arr6[j] = y6;

         y7 = arr7[i];
         arr7[i] = arr7[j];
         arr7[j] = y7;

         y8 = arr8[i];
         arr8[i] = arr8[j];
         arr8[j] = y8;

         i++;
         j--;
      }
   } while (i <= j);

   // recurse
   if (low < j)
      quicksort(arr0, arr1, arr2, arr3, arr4, arr5, arr6, arr7, arr8, low, j);

   if (i < high)
      quicksort(arr0, arr1, arr2, arr3, arr4, arr5, arr6, arr7, arr8, i, high);
};

ClassImp(AtSOLARISUnpacker);