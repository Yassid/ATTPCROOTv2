#include "AtMAGNEXParsingTask.h"

#include "AtEvent.h"
#include "AtHit.h"

#include <FairLogger.h>      // for LOG
#include <FairRootManager.h> // for FairRootManager

#include <TObject.h> // for TObject

#include <H5Cpp.h>

#include <utility> // for move

AtMAGNEXParsingTask::AtMAGNEXParsingTask(TString fileName, TString mapFileName, TString outputBranchName)
   : fInputFileName(std::move(fileName)), fInputMapFileName(std::move(mapFileName)),
     fOutputBranchName(std::move(outputBranchName)), fEventArray("AtEvent", 1)
{
}

InitStatus AtMAGNEXParsingTask::Init()
{

   FairRootManager *ioMan = FairRootManager::Instance();
   if (ioMan == nullptr) {
      LOG(ERROR) << "Cannot find RootManager!";
      return kERROR;
   }

   ioMan->Register(fOutputBranchName, "AtTPC", &fEventArray, fIsPersistence);

   fFile = std::make_unique<TFile>(fInputFileName, "READ");
   if (fFile->IsZombie()) {
      LOG(ERROR) << "Cannot open ROOT file " << fInputFileName;
      return kERROR;
   } else {
      std::cout << "File " << fInputFileName << " opened successfully"
                << "\n";
   }

   fMapFile = std::make_unique<std::ifstream>(fInputMapFileName);
   fPadMap = std::make_unique<std::map<Int_t, Int_t>>();

   if (!fMapFile->is_open()) {
      LOG(ERROR) << "Cannot open map file " << fInputMapFileName;
      return kERROR;
   } else {
      std::cout << "Map file " << fInputMapFileName << " opened successfully"
                << "\n";

      std::string line;
      for (auto i = 0; i < 8; ++i) {
         std::getline(*fMapFile, line, '\n');
         LOG(debug) << line << "\n";
      }

      while (!fMapFile->eof()) {
         Int_t channel, pad;
         std::getline(*fMapFile, line, '\n');
         std::istringstream buffer(line);
         buffer >> channel >> pad;
         LOG(debug) << "Channel: " << channel << " Pad: " << pad << "\n";
         fPadMap->insert(std::make_pair(channel, pad));
      }
   }

   fTree = (TTree *)fFile->Get("Data_R");
   std::cout << "Number of entries in the tree: " << fTree->GetEntries() << "\n";
   fTree->SetBranchAddress("Board", &Board);
   fTree->SetBranchAddress("Channel", &Channel);
   fTree->SetBranchAddress("FineTSInt", &FTS);
   fTree->SetBranchAddress("CoarseTSInt", &CTS);
   fTree->SetBranchAddress("Timestamp", &Timestamp);
   fTree->SetBranchAddress("Charge", &Charge);
   fTree->SetBranchAddress("Flags", &Flags);
   fTree->SetBranchAddress("Pads", &pad);
   fTree->SetBranchAddress("Charge_cal", &Charge_cal);
   fTree->SetBranchAddress("Row", &Row);
   fTree->SetBranchAddress("Section", &Section);

   return kSUCCESS;
}

void AtMAGNEXParsingTask::Exec(Option_t *opt)
{

   fTree->GetEntry(fEntryNum);
   ++fEntryNum;
   ULong64_t timeinit = Timestamp;
   std::cout << " time init: " << timeinit << "\n";
   ULong64_t hitNum = 0;

   auto *event = dynamic_cast<AtEvent *>(fEventArray.ConstructedAt(0, "C")); // Get and clear old event
   event->SetEventID(fEventNum++);
   event->SetTimestamp(Timestamp);

   // First hit
   auto pad = fPadMap->find(Channel)->second;
   Double_t x = 2.5 + 5 * pad;
   Double_t y = FTS; // Needs calibration
   Double_t z = Row * 21.2 + 18.60;
   auto &fhit = event->AddHit(hitNum, AtHit::XYZPoint(x, y, z), Charge_cal);
   fhit.SetHitID(hitNum);
   fhit.SetPadNum(pad);
   fhit.SetTimeStamp(FTS);
   ++hitNum;

   // std::cout<<fhit.GetHitID()<<" "<<fhit.GetPadNum()<<" "<<fhit.GetTimeStamp()<<" "<<fhit.GetCharge()<<"
   // "<<fEntryNum<<"\n";

   while (fEntryNum < fTree->GetEntries()) {

      fTree->GetEntry(fEntryNum);
      ++fEntryNum;
      // std::cout<<Timestamp<<"   "<<Timestamp - timeinit<<"\n";
      if (Timestamp - timeinit > fWindowSize) {
         --fEntryNum;
         break;
      }

      pad = fPadMap->find(Channel)->second;
      x = 2.5 + 5 * pad;
      y = FTS; // Needs calibration
      z = Row * 21.2 + 18.60;
      auto &hit = event->AddHit(hitNum, AtHit::XYZPoint(x, 0, z), Charge_cal);
      hit.SetHitID(hitNum);
      hit.SetPadNum(pad);
      hit.SetTimeStamp(FTS);
      ++hitNum;

      // std::cout<<hit.GetHitID()<<" "<<hit.GetPadNum()<<" "<<hit.GetTimeStamp()<<" "<<hit.GetCharge()<<"
      // "<<fEntryNum<<"\n";
   }

   // std::cout<<" Out of loop "<<Timestamp<<"   "<<Timestamp - timeinit<<"\n";

   std::cout << " Event number : " << fEventNum << " Entry number " << fEntryNum << " Number of hits : " << hitNum
             << "\n";
}

ClassImp(AtMAGNEXParsingTask);
