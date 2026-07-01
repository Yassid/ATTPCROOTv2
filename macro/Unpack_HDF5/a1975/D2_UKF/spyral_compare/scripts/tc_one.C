R__LOAD_LIBRARY(libAtData)
R__LOAD_LIBRARY(libAtReconstruction)
#include "AtTrackFinderTC.h"
#include "AtEvent.h"
#include "AtPatternEvent.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
#include "TStopwatch.h"
void tc_one(int ev=49, float s=0.3,int k=19,int n=2,int m=15,float r=2,float a=0.03,float t=4.0){
  TFile* f=TFile::Open("/mnt/f/a1975/reco_d2/run_0016_reco.root");
  TTree* tr=(TTree*)f->Get("cbmsim");
  TClonesArray* corr=nullptr; tr->SetBranchAddress("AtEventCorrected",&corr);
  tr->GetEntry(ev); AtEvent* e=(AtEvent*)corr->At(0);
  TStopwatch sw; sw.Start();
  AtPATTERN::AtTrackFinderTC tf;
  tf.SetScluster(s);tf.SetKtriplet(k);tf.SetNtriplet(n);tf.SetMcluster(m);tf.SetRsmooth(r);tf.SetAtriplet(a);tf.SetTcluster(t);
  auto pe=tf.FindTracks(*e);
  int nh=e->GetNumHits(),ntr=0,assigned=0,maxt=0;
  for(auto& x:pe->GetTrackCand()){int hn=x.GetHitArray().size();assigned+=hn;maxt=std::max(maxt,hn);ntr++;}
  printf("RESULT ev%d s=%.2f k=%d m=%d a=%.2f t=%.1f : ntracks=%d maxtrack=%d assigned=%d/%d (%.0f%%) [%.1fs]\n",
    ev,s,k,m,a,t,ntr,maxt,assigned,nh,100.*assigned/nh,sw.RealTime());
}
