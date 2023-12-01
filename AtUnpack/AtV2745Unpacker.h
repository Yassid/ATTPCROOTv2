/* V2745 CAEN Digitizer unpacker for MAGNEX tracker
 *
 *  Author: Yassid Ayyad (yassid.ayyad@usc.es)
 *  Log: Class started on 19 July 2023
 */

#ifndef _ATV2745UNPACKER_H_
#define _ATV2745UNPACKER_H_

#include "AtUnpacker.h"

#include <Rtypes.h>
#include <TString.h>

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

using dataidpair = std::pair<std::ifstream, int>;
using datanID = std::vector<std::unique_ptr<dataidpair>>;

class AtV2745Unpacker : public AtUnpacker {
protected:
   std::mutex fRawEventMutex;
   Int_t fNumFiles;

   AtPad *createPad(const struct data &block);

public:
   AtV2745Unpacker(mapPtr map, Long64_t numBlocks = 1E3, Int_t numFiles = 5);
   ~AtV2745Unpacker() = default;

   virtual void Init() override;
   virtual void FillRawEvent(AtRawEvent &event) override; // Pass by ref to ensure it's a valid object
   virtual bool IsLastEvent() override;
   virtual void SetInputFileName(std::string fileName) override;

   virtual Long64_t GetNumEvents() override;

   void ProcessFile(dataidpair *datawithid);
   void PreprocessData();

private:
   std::string fInputFileName{};
   datanID fBinaryDataFiles;
   Long64_t fNumBlocks; // Number of blocks from V2745 to unpack

   const Double_t us_to_ps = 1.0e6;               // Conversion factor micros -> ps
   const Double_t fineTS_to_ps = 8000.0 / 1024.0; // Conversion factor fine timestamp -> ps (8 ns over 1024 bins)

   void ProcessInputFile();
   Int_t GetBoardID(std::string fileName);

   ClassDefOverride(AtV2745Unpacker, 1)
};

/*
 * Struct used as single channel data record in binary files. Has to be identical to struct used in DAQ program
 */
struct data {
   uint8_t Channel;
   Double_t Coarse_Time_micros;
   uint16_t Fine_Time_int;
   uint16_t LoPFlags;
   uint16_t HiPFlags;
   uint16_t Energy;
};

#endif //#ifndef _ATV2745UNPACKER_H_
