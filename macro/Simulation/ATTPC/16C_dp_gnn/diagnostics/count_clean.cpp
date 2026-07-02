#include "AtEvent.h"
#include "AtHit.h"
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <cstdio>
int main(int argc, char** argv){
  TFile* f=TFile::Open(argv[1]); TTree* t=(TTree*)f->Get("cbmsim");
  TClonesArray *H=nullptr,*C=nullptr;
  t->SetBranchAddress("AtEventH",&H); t->SetBranchAddress("AtEventClean",&C);
  long th=0,tc=0,nev=0;
  for(long i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(H&&H->GetEntriesFast()){auto*e=(AtEvent*)H->At(0); th+=e->GetHits().size();}
    if(C&&C->GetEntriesFast()){auto*e=(AtEvent*)C->At(0); tc+=e->GetHits().size();}
    nev++;
  }
  printf("events %ld | AtEventH hits %ld (%.0f/ev) | AtEventClean hits %ld (%.0f/ev) | kept %.1f%%, removed %.1f%%\n",
         nev, th, (double)th/nev, tc, (double)tc/nev, 100.0*tc/th, 100.0*(th-tc)/th);
  return 0;
}
