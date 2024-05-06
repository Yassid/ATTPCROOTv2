#ifndef ATMAGNEXPARSINGTASK_H
#define ATMAGNEXPARSINGTASK_H

#include <FairTask.h>

#include <Rtypes.h> // for Bool_t, Int_t, Option_t
#include <TClonesArray.h>
#include <TFile.h>   // for TFile
#include <TString.h> // for TString

#include <H5Cpp.h>

#include <memory> // for unique_ptr

class TBuffer;
class TClass;
class TMemberInspector;

class AtMAGNEXParsingTask : public FairTask {

protected:
   Bool_t fIsPersistence{false};
   TString fInputFileName;
   TString fInputMapFileName;
   TString fOutputBranchName;

   std::unique_ptr<TFile> fFile{nullptr}; //!
   TTree *fTree{nullptr};                 //!
   TClonesArray fEventArray;
   std::unique_ptr<std::ifstream> fMapFile;
   std::unique_ptr<std::map<Int_t, Int_t>> fPadMap;

   Int_t fEventNum{0};
   ULong64_t fEntryNum{0};
   ULong64_t fWindowSize{2000000};

   UShort_t Channel;
   UShort_t pad;
   UShort_t FTS;
   ULong64_t CTS;
   ULong64_t Timestamp;
   UShort_t Board;
   UShort_t Charge;
   Double_t Charge_cal;
   UInt_t Flags;
   UShort_t Row;
   UShort_t Section;

public:
   AtMAGNEXParsingTask(TString fileName, TString mapFileName, TString outputBranchName = "AtEventH");

   void SetPersistence(bool val) { fIsPersistence = val; }

   virtual InitStatus Init() override;
   virtual void Exec(Option_t *opt) override;

   ClassDefOverride(AtMAGNEXParsingTask, 1);
};

#endif //#ifndef ATMAGNEXPARSINGTASK_H
