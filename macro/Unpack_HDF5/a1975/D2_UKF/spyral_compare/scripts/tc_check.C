R__LOAD_LIBRARY(libAtData)
R__LOAD_LIBRARY(libAtReconstruction)
#include "AtTrackFinderTC.h"
#include "AtEvent.h"
#include "AtPatternEvent.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
static void run(AtEvent* e,float a,float t,const char* tag){
  AtPATTERN::AtTrackFinderTC tf;
  tf.SetClusterRadius(15);tf.SetClusterDistance(7.5);tf.SetUseSelectAndMerge(true);
  tf.SetScluster(0.3);tf.SetKtriplet(19);tf.SetNtriplet(2);tf.SetMcluster(15);tf.SetRsmooth(2);tf.SetAtriplet(a);tf.SetTcluster(t);
  auto pe=tf.FindTracks(*e);
  int nh=e->GetNumHits(),ntr=0,asg=0,maxt=0;
  for(auto& x:pe->GetTrackCand()){int hn=x.GetHitArray().size();asg+=hn;maxt=std::max(maxt,hn);ntr++;}
  printf(">> %-14s ntr=%2d max=%4d asg=%4d/%d (%3.0f%%)\n",tag,ntr,maxt,asg,nh,100.*asg/nh);
}
void tc_check(int ev){
  TFile* f=TFile::Open("/mnt/f/a1975/reco_d2/run_0016_reco.root");
  TTree* tr=(TTree*)f->Get("cbmsim"); TClonesArray* corr=nullptr;
  tr->SetBranchAddress("AtEventCorrected",&corr); tr->GetEntry(ev); AtEvent* e=(AtEvent*)corr->At(0);
  printf("== ev %d (%d hits) ==\n",ev,e->GetNumHits());
  run(e,0.03,4.0,"DEFAULT"); run(e,0.20,12.0,"a=.20 t=12");
}
