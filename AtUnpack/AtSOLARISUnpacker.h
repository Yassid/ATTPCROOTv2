#ifndef _ATSOLARISUNPACKER_H_
#define _ATSOLARISUNPACKER_H_

#include "AtUnpacker.h"

#include <Rtypes.h>
#include <TString.h>

#include "SolReader.h"

#include <memory>
#include <mutex>
#include <regex> // for match_results<>::_Base_type, regex_replace
#include <string>
#include <utility> // for pair
#include <vector>

class AtRawEvent;
class AtPad;
class TBuffer;
class TClass;
class TMemberInspector;

class AtSOLARISUnpacker : public AtUnpacker {
protected:
   Int_t fNumFiles;

   int verbosity = 0;
   char fileCalCharge[200] = "Calib_files/chargeCalib_PA.txt"; // charge calibration file
   char filePadMapping[200] = "Calib_files/channel2pad_3.txt"; // pads mapping file
   Int_t board_id[20] = {22642, 22643, 22644, 22645, 21247};   // board id for each file

   Double_t cal_fac[5][64]; // matrix 5x64 for calibration factors for 5 PAs x 64 channels
   UShort_t pad[64];        // the i-th value of the array is the pad corresponding to the the i-th

public:
   AtSOLARISUnpacker(mapPtr map, Long64_t numBlocks = 1E3, std::string fileName = "");
   ~AtSOLARISUnpacker() = default;

   virtual void Init() override;
   virtual void FillRawEvent(AtRawEvent &event) override; // Pass by ref to ensure it's a valid object
   virtual bool IsLastEvent() override;
   virtual void SetInputFileName(std::string fileName) override;
   virtual Long64_t GetNumEvents() override;

   void ProcessFile();
   void SetWindowSize(ULong64_t windowSize);

private:
   std::string fInputFileName{};
   Long64_t fNumBlocks;            // Number of blocks from V2745 to unpack
   ULong64_t fWindowSize{2000000}; // in ps

   const Double_t us_to_ps = 1.0e6;               // Conversion factor micros -> ps
   const Double_t ns_to_ps = 1.0e+03;             // Conversion factor nanos -> ps
   const Double_t fineTS_to_ps = 8000.0 / 1024.0; // Conversion factor fine timestamp -> ps (8 ns over 1024 bins)

   void ProcessInputFile();
   void processData();

   Int_t findRow(int dig, int pad);
   Int_t findColumn(int dig);
   Int_t findSection(int dig);
   void calib_init(int ver = 0);
   void pad_init(int ver = 0);
   void quicksort(UShort_t arr0[], ULong64_t arr1[], UShort_t arr2[], UShort_t arr3[], UShort_t arr4[],
                  ULong64_t arr5[], UInt_t arr6[], Double_t arr7[], UShort_t arr8[], int low, int high);

   std::unique_ptr<SolReader> fReader;
   std::unique_ptr<Hit> hit;

   ClassDefOverride(AtSOLARISUnpacker, 1)
};

#endif //#ifndef _ATSOLARISUNPACKER_H_
