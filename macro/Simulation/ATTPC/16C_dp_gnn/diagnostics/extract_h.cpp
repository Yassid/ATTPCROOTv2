#include "AtEvent.h"
#include "AtHit.h"
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <cstdio>
int main(int argc,char**argv){
  TFile*f=TFile::Open(argv[1]); TTree*t=(TTree*)f->Get("cbmsim");
  TClonesArray*H=nullptr; t->SetBranchAddress(argv[3],&H);
  FILE*fp=fopen(argv[2],"w"); fprintf(fp,"event,x,y,z,q\n");
  for(long i=0;i<t->GetEntries();++i){t->GetEntry(i);
    if(!H||!H->GetEntriesFast())continue; auto*e=(AtEvent*)H->At(0);
    for(auto&h:e->GetHits()){auto p=h->GetPosition();
      fprintf(fp,"%ld,%.2f,%.2f,%.2f,%.2f\n",i,p.X(),p.Y(),p.Z(),h->GetCharge());}}
  fclose(fp); return 0;
}
