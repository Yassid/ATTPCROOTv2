#include "AtMAGNEXParserTask.h"

#include "AtEvent.h"
#include "AtHit.h"

#include <FairLogger.h>      // for LOG
#include <FairRootManager.h> // for FairRootManager

#include <TObject.h> // for TObject

#include <utility> // for move

AtMAGNEXParserTask::AtMAGNEXParserTask(TString fileName, TString branchName)
   : fInputFileName(std::move(fileName)), fOutputBranchName(std::move(branchName))
{
}

InitStatus AtMAGNEXParserTask::Init()
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
   }

   fTree = (TTree *)fFile->Get("Data_R");

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

void AtMAGNEXParserTask::Exec(Option_t *opt)
{
   fTree->GetEntry(0);
   ULong64_t timeinit = Timestamp;
   std::cout << " time init: " << timeinit << "\n";

   auto *event = dynamic_cast<AtEvent *>(fEventArray.ConstructedAt(0, "C")); // Get and clear old event

   AtHit_t hitFromFile{};
   event->AddHit(-1, AtHit::XYZPoint(hitFromFile.x, hitFromFile.y, hitFromFile.z), hitFromFile.A);

   // After this loop filling the event, I don't think there is anything else to do
}

ClassImp(AtMAGNEXParserTask);
