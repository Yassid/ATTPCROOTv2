#include "AtTrackingEvent.h"
#include "AtFittedTrack.h"
#include <TClonesArray.h>
#include <TFile.h>
#include <TTree.h>
#include <cstdio>
int main(int argc,char**argv){
  TFile*f=TFile::Open(argv[1]); TTree*t=(TTree*)f->Get("cbmsim");
  TClonesArray*trk=nullptr; t->SetBranchAddress("AtTrackingEvent",&trk);
  FILE*fp=fopen(argv[2],"w"); fprintf(fp,"event,ke,theta,chi2ndf\n");
  long nrow=0;
  for(long i=0;i<t->GetEntries();++i){ t->GetEntry(i);
    if(!trk||!trk->GetEntriesFast())continue; auto*ev=(AtTrackingEvent*)trk->At(0); if(!ev)continue;
    for(auto&ft:ev->GetFittedTracks()){ if(!ft||!ft->GetTrackMetadata())continue;
      double ke,th,ndf,chi2;
      try { auto&k=ft->GetKinematics(); ke=k.kineticEnergy; th=k.theta;
            ndf=ft->GetTrackMetadata()->GetNdf(); chi2=ft->GetTrackMetadata()->GetChi2(); }
      catch(...){ continue; }
      double c2n=ndf>0?chi2/ndf:1e9;
      fprintf(fp,"%ld,%.4f,%.5f,%.3f\n",i,ke,th,c2n); nrow++; } }
  fclose(fp); fprintf(stderr,"wrote %ld rows\n",nrow); return 0;
}
