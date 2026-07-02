#include "AtTrackingEvent.h"
#include "AtFittedTrack.h"
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <cstdio>
int main(int argc,char**argv){
  TFile*f=TFile::Open(argv[1]); TTree*t=(TTree*)f->Get("cbmsim");
  TClonesArray*trk=nullptr; t->SetBranchAddress("AtTrackingEvent",&trk);
  FILE*fp=fopen(argv[2],"w"); fprintf(fp,"event,trackid,chi2,ndf\n");
  for(long i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!trk||!trk->GetEntriesFast())continue; auto*ev=(AtTrackingEvent*)trk->At(0); if(!ev)continue;
    for(auto&ft:ev->GetFittedTracks()){ if(!ft||!ft->GetTrackMetadata())continue;
      double ndf=ft->GetTrackMetadata()->GetNdf(), chi2=ft->GetTrackMetadata()->GetChi2();
      fprintf(fp,"%ld,%d,%.3f,%.1f\n",i,ft->GetTrackID(),chi2,ndf);} }
  fclose(fp); return 0;
}
