R__LOAD_LIBRARY(libAtData)
R__LOAD_LIBRARY(libAtReconstruction)
#include "AtTrackFinderTC.h"
#include "AtEvent.h"
#include "AtPatternEvent.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"
static AtEvent* gE=nullptr;
static void cfg(float cr,float cd,bool sel,float s,int k,int m,float a,float t,const char* tag){
  AtPATTERN::AtTrackFinderTC tf;
  tf.SetClusterRadius(cr); tf.SetClusterDistance(cd); tf.SetUseSelectAndMerge(sel);
  tf.SetScluster(s);tf.SetKtriplet(k);tf.SetNtriplet(2);tf.SetMcluster(m);tf.SetRsmooth(2);tf.SetAtriplet(a);tf.SetTcluster(t);
  auto pe=tf.FindTracks(*gE);
  int nh=gE->GetNumHits(),ntr=0,assigned=0,maxt=0;
  for(auto& x:pe->GetTrackCand()){int hn=x.GetHitArray().size();assigned+=hn;maxt=std::max(maxt,hn);ntr++;}
  printf(">> %-26s ntr=%2d max=%4d asg=%4d/%d (%3.0f%%)\n",tag,ntr,maxt,assigned,nh,100.*assigned/nh);
}
void tc_sweep(int ev=49){
  TFile* f=TFile::Open("/mnt/f/a1975/reco_d2/run_0016_reco.root");
  TTree* tr=(TTree*)f->Get("cbmsim"); TClonesArray* corr=nullptr;
  tr->SetBranchAddress("AtEventCorrected",&corr); tr->GetEntry(ev); gE=(AtEvent*)corr->At(0);
  printf("==== event %d : cloud %d hits ====\n",ev,gE->GetNumHits());
  cfg(15,7.5,true, 0.3,19,15,0.03,4.0, "DEFAULT(reco)");
  cfg(15,7.5,true, 0.3,19,15,0.10,8.0, "a=.10 t=8");
  cfg(15,7.5,true, 0.3,19,15,0.20,12.0,"a=.20 t=12");
  cfg(30,15, true, 0.3,19,15,0.10,8.0, "merge:cr30 cd15 a.10 t8");
  cfg(40,20, true, 0.3,19,8, 0.20,12.0,"loose:cr40 cd20 m8 a.20 t12");
  cfg(15,7.5,false,0.3,19,15,0.03,4.0, "no-select default");
  cfg(15,7.5,false,0.3,19,8, 0.20,12.0,"no-select loose");
}
