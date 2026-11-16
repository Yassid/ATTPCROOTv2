#include "resfitC.C"
// The Ebeam(vz) calibration is flat above vz ~ 300 and breaks below it, in BOTH codes: the g.s.
// centroid walks by 0.21 MeV (ATTPCROOT) / 0.12 (Spyral) and the peak broadens by 10-30 %.
// Spyral's adopted window was 0-400, i.e. mostly the BAD region -- which is why the two codes
// disagreed about the continuum.  Adopt vz in [300,1000] for both and see what follows.
void vzadopt(){
  setCM(10,5,60);
  RH::exLo=3.0; RH::exHi=9.0; RH::TIE12=false; RH::ONETERM=true;
  for(int i=0;i<RH::NR;i++) RH::USE[i]=true;
  Cfg b0=CFG[0], b1=CFG[1];
  CFG[0].vzHi=1000; CFG[1].vzHi=1000; VZLO=300; BB_QUIET=true; binbuild();
  for(int ic=0;ic<2;ic++){
    TFile f("dd_cmbins.root"); double N=0; TH1D* hb=nullptr; double bN=-1;
    for(int b=0;b<NCM;b++){ TH1D* x=(TH1D*)f.Get(Form("%s_cm%d",CFG[ic].nm,b)); if(!x) continue;
      N+=x->GetEntries(); if(x->GetEntries()>bN){bN=x->GetEntries(); hb=(TH1D*)x->Clone("hb"); hb->SetDirectory(nullptr);} }
    f.Close();
    double mu=0,fw=0;
    if(hb){ TF1 g(Form("gq%d",ic),"gaus",-1.0,1.0);
            g.SetParameters(hb->GetBinContent(hb->GetMaximumBin()),0,0.25);
            hb->Fit(&g,"RQ0N"); mu=g.GetParameter(1); fw=2.355*g.GetParameter(2); }
    RH::BKG=RH::HALON; double vf=rhRun(CFG[ic],NCM,false), nh=RH::LASTP0;
    double s=0; { RH::SBFIX=true;
      for(double d=0.02;d<1.5;d+=0.02){ RH::SBVAL=nh+d; if(rhRun(CFG[ic],NCM,false)-vf>=1.0){s=d;break;} }
      RH::SBFIX=false; }
    RH::BKG=RH::HALO; double vH=rhRun(CFG[ic],NCM,false);
    RH::BKG=RH::EXP;  double vE=rhRun(CFG[ic],NCM,false);
    RH::BKG=RH::PS;   double vP=rhRun(CFG[ic],NCM,false);
    printf("\n================ %s   vz in [300,1000] ================\n",CFG[ic].nm);
    printf("   N = %.0f     g.s. centroid %.3f   FWHM %.3f\n",N,mu,fw);
    printf("   halo power n = %.3f +- %.3f      (%.1f sigma from the predicted 2)\n",nh,s,
           std::fabs(nh-2.0)/std::max(1e-6,s));
    printf("   -2lnL   halo n=2 (0 par) %8.1f | free n (1 par) %8.1f | exp (1 par) %8.1f | uniform PS (0 par) %8.1f\n",
           vH,vf,vE,vP);
    printf("   parameter-free halo vs empirical exponential: %+.1f with ONE FEWER parameter\n",vH-vE);
  }
  CFG[0]=b0; CFG[1]=b1; VZLO=0; BB_QUIET=true; binbuild();
}
