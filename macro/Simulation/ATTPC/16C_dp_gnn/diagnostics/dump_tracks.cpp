#include "AtPatternEvent.h"
#include "AtTrack.h"
#include "AtHit.h"
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <cstdio>
int main(int argc,char**argv){
  TFile*f=TFile::Open(argv[1]); TTree*t=(TTree*)f->Get("cbmsim");
  TClonesArray*arr=nullptr; t->SetBranchAddress("AtPatternEvent",&arr);
  FILE*fp=fopen(argv[2],"w"); fprintf(fp,"event,trackid,x,y,z,q\n");
  for(long i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!arr||!arr->GetEntriesFast())continue; auto*pe=(AtPatternEvent*)arr->At(0); if(!pe)continue;
    auto&tracks=pe->GetTrackCand(); int tid=0;
    for(auto&tr:tracks){ for(auto&h:tr.GetHitArray()){ if(!h)continue; auto p=h->GetPosition();
      fprintf(fp,"%ld,%d,%.2f,%.2f,%.2f,%.2f\n",i,tid,p.X(),p.Y(),p.Z(),h->GetCharge()); } tid++; } }
  fclose(fp); return 0;
}
