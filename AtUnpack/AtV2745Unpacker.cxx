#include "AtV2745Unpacker.h"

#include "AtMap.h"
#include "AtPad.h"
#include "AtPadBase.h" // for AtPadBase
#include "AtPadV2745.h"
#include "AtPadValue.h"
#include "AtRawEvent.h"
#include "AtTpcMap.h"

#include <FairLogger.h>

#include <Rtypes.h>

#include <algorithm>
#include <array> // for array
#include <future>
#include <iostream>
#include <iterator> // for begin, end
#include <numeric>  // for accumulate
#include <thread>
#include <utility>

AtV2745Unpacker::AtV2745Unpacker(mapPtr map, Long64_t numBlocks, Int_t numFiles) : AtUnpacker(map), fNumFiles(numFiles)
{
   fNumBlocks = numBlocks;
}

void AtV2745Unpacker::Init()
{
   fNumFiles = 5;
   ProcessInputFile();

   // Unpack blocks for event building
   for (auto i = 0; i < fNumBlocks; ++i)
      PreprocessData();
}

void AtV2745Unpacker::PreprocessData()
{
   std::vector<std::thread> fileProcessors;

   for (Int_t iFile = 0; iFile < fNumFiles; iFile++) {
      fileProcessors.emplace_back([this](Int_t fileIdx) { this->ProcessFile(fBinaryDataFiles[fileIdx].get()); }, iFile);
   }

   for (auto &fileThread : fileProcessors)
      fileThread.join();
}

void AtV2745Unpacker::FillRawEvent(AtRawEvent &event)
{
   fRawEvent = &event;
}

void AtV2745Unpacker::ProcessFile(dataidpair *datawithid)
{
   auto datastr = &datawithid->first;
   auto dataid = datawithid->second;

   LOG(debug) << " Processing data from board " << dataid << "\n";

   char block[sizeof(struct data)]; // Block to read data from bin file (requires char* type)

   struct data *data; // Single channel data from file

   // Read entry
   datastr->read(block, sizeof(struct data));
   data = (struct data *)block;

   std::lock_guard<std::mutex> lk(fRawEventMutex);
   LOG(info) << dataid << "  " << (UShort_t)(data->Channel) << " " << (ULong64_t)(data->Coarse_Time_micros * us_to_ps)
             << "  " << (UShort_t)(data->Energy) << "   "
             << "i"
             << "\n";
}

AtPad *AtV2745Unpacker::createPad(const struct data &block)
{

   AtV2745Block v2745_block;
   v2745_block.channel = (UShort_t)(block.Channel);
   v2745_block.charge = (UShort_t)(block.Energy);
   v2745_block.time_ps =
      (ULong64_t)(us_to_ps * block.Coarse_Time_micros) + ULong64_t((Double_t)(block.Fine_Time_int) * fineTS_to_ps);
   v2745_block.coarse_time_int = block.Coarse_Time_micros;
   v2745_block.fine_time_int = (UShort_t)(block.Fine_Time_int);
   auto LoPFlags0 = (UInt_t)(block.LoPFlags);
   auto HiPFlags0 = (UInt_t)(block.HiPFlags);
   v2745_block.flags = LoPFlags0 << 16;
   v2745_block.flags += HiPFlags0;

   // TODO BoardID, row, column

   auto pad = fRawEvent->AddAuxPad(std::to_string(v2745_block.channel)).first;
   pad->AddAugment("Block", std::make_unique<AtPadV2745>(v2745_block));
   return pad;
}

bool AtV2745Unpacker::IsLastEvent()
{
   return 0;
}

void AtV2745Unpacker::SetInputFileName(std::string fileName)
{
   fInputFileName = fileName;
}

Long64_t AtV2745Unpacker::GetNumEvents()
{

   return 0;
}

void AtV2745Unpacker::ProcessInputFile()
{
   LOG(info) << " Processing input file : " << fInputFileName << "\n";

   std::ifstream listFile(fInputFileName);
   TString dataFileWithPath;

   Int_t numFiles;

   while (dataFileWithPath.ReadLine(listFile)) {
      auto id = GetBoardID(dataFileWithPath.Data());

      LOG(info) << " Adding file : " << dataFileWithPath << " for board ID " << id << "\n";

      std::ifstream dataFile(dataFileWithPath.Data(), std::ios::binary);

      if (dataFile.fail()) {
         std::cerr << " File not found. Exiting... "
                   << "\n";
         std::exit(0);
      }

      std::unique_ptr<dataidpair> datanid_ptr = std::make_unique<dataidpair>(std::make_pair(std::move(dataFile), id));
      fBinaryDataFiles.push_back(std::move(datanid_ptr));
      ++numFiles;
   }
}

Int_t AtV2745Unpacker::GetBoardID(std::string fileName)
{
   std::string input_string = fileName;

   std::regex number_pattern("\\d{5}");

   std::sregex_iterator iterator(input_string.begin(), input_string.end(), number_pattern);
   std::sregex_iterator end_iterator;
   std::smatch match;

   while (iterator != end_iterator) {
      match = *iterator;
      LOG(debug) << "Found Board ID: " << match.str() << std::endl;
      ++iterator;
   }

   return std::stoi(match.str());
}

ClassImp(AtV2745Unpacker);
