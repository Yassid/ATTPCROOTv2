#include "binacc.h"
#include "binbuild.C"
// Loosening chi2 buys statistics -- but the state list depends on RESOLUTION (two pairs are
// separated by 0.50 and 0.70 against FWHM 0.93).  If a looser cut degrades the g.s. width the
// trade is a bad one.  Scan it and watch both.
void chi2opt(){
  for(int ic=0;ic<2;ic++){
    Cfg base=CFG[ic];
    double cuts[9]={0.1,0.2,0.3,0.5,1.0,2.0,5.0,20.0,1e9};
    printf("\n================ %s : chi2 cut vs statistics AND resolution ================\n",base.nm);
    printf("   %-10s %10s %8s %10s %10s %10s\n","chi2 <","N(5-9)","gain","gs FWHM","gs centr","peak/val");
    double N0=-1;
    for(double cut:cuts){
      Cfg k=base; k.chi2=cut; k.vzHi=900;
      Pars A=parsOf(k); A.vzLo=300;
      TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
      double ke,th,vz,c2,icd=0; float icf=0;
      t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
      t->SetBranchAddress("chi2ndf",&c2);
      bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
      if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
      TH1D* hg=new TH1D(Form("g%d%.0f",ic,cut*10),"",160,-2,2); hg->SetDirectory(nullptr);
      long n59=0,npk=0,nvl=0;
      for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
        double icv=icF?icf:icd;
        if(c2>cut||ke<k.keLo||ke>k.keHi||vz<300||vz>900||icv<k.icLo||icv>k.icHi) continue;
        double thf=std::fabs(th); if(thf>90) thf=180-thf;
        double eb=ebAt(A,vz,true), thc=thCorr(A,thf,ke);
        double e=exOf(A,eb,thc,ke); if(!std::isfinite(e)) continue;
        double tcm=thcmProj(A,eb,thc,ke); if(tcm<10||tcm>60) continue;
        hg->Fill(e);
        if(e>5&&e<9) n59++;
        if(e>6.4&&e<7.0) npk++;
        if(e>7.8&&e<8.2) nvl++; }
      f->Close();
      TF1 g(Form("gg%d%.0f",ic,cut*10),"gaus",-1.0,1.0);
      g.SetParameters(hg->GetBinContent(hg->GetMaximumBin()),0,0.25);
      hg->Fit(&g,"RQ0N");
      if(N0<0) N0=n59;
      printf("   %-10.4g %10ld %7.2fx %10.3f %10.3f %10.3f\n",
             cut,n59,N0>0?(double)n59/N0:0,2.355*g.GetParameter(2),g.GetParameter(1),
             nvl>0?(double)npk/nvl*(0.4/0.6):0.0); }
  }
}
