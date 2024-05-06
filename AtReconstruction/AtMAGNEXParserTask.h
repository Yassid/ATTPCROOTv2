#ifndef ATMAGNEXPARSERTASK_H
#define ATMAGNEXPARSERTASK_H

#include <FairTask.h>

#include <Rtypes.h> // for Bool_t, Int_t, Option_t
#include <TClonesArray.h>
#include <TFile.h>   // for TFile
#include <TString.h> // for TString
#include <TTree.h>   // for TTree

#include <memory> // for unique_ptr

class TBuffer;
class TClass;
class TMemberInspector;

class AtMAGNEXParserTask : public FairTask {

protected:
   Bool_t fIsPersistence{false};
   TString fInputFileName;
   TString fOutputBranchName;

   std::unique_ptr<TFile> fFile{nullptr}; //!
   TTree *fTree{nullptr};                 //!
   TClonesArray fEventArray;

   Int_t fEventNum{0};

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
   AtMAGNEXParserTask(TString fileName, TString outputBranchName = "AtEventH");

   void SetPersistence(bool val) { fIsPersistence = val; }

   virtual InitStatus Init() override;
   virtual void Exec(Option_t *opt) override;

   ClassDefOverride(AtMAGNEXParserTask, 1);
};

#endif //#ifndef ATMAGNEXPARSERTASK_H