#include "binacc.h"
#include "binbuild.C"
// Decisive calibration test before opening Spyral's vz window: is the g.s. centroid FLAT against
// vz all the way out to 1000?  If Ebeam(vz) had the wrong sign or slope the centroid would walk.
void ebzcheck(){
  for(int ic=0;ic<2;ic++){
    Cfg& k=CFG[ic]; Pars A=parsOf(k); A.vzLo=-1e9; A.vzHi=1e9;
    TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
    double ke,th,vz,c2,icd=0; float icf=0;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
    t->SetBranchAddress("chi2ndf",&c2);
    bool icF = TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
    if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
    const int NS=10; double lo=0,hi=1000; TH1D* h[NS];
    for(int s=0;s<NS;s++){ h[s]=new TH1D(Form("h%d%d",ic,s),"",160,-2,2); h[s]->SetDirectory(nullptr); }
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      double icv = icF?icf:icd;
      if(c2>k.chi2||ke<k.keLo||ke>k.keHi||icv<k.icLo||icv>k.icHi) continue;
      if(vz<lo||vz>=hi) continue;
      double thf=std::fabs(th); if(thf>90) thf=180-thf;
      double thc=thCorr(A,thf,ke), eb=ebAt(A,vz,true);
      double e=exOf(A,eb,thc,ke); if(!std::isfinite(e)) continue;
      int idx=(int)((vz-lo)/((hi-lo)/NS)); if(idx<0||idx>=NS) continue;
      h[idx]->Fill(e); }
    f->Close();
    printf("\n================ %s : g.s. centroid vs vz  (Ebeam(vz) ON) ================\n",k.nm);
    printf("   %-12s %8s %10s %9s\n","vz slice","N","centroid","FWHM");
    double mn=1e9,mx=-1e9;
    for(int s=0;s<NS;s++){ if(h[s]->GetEntries()<200){ printf("   %4.0f-%4.0f    %8.0f      (thin)\n",
         lo+s*100.0,lo+(s+1)*100.0,h[s]->GetEntries()); continue; }
      TF1 g(Form("gg%d%d",ic,s),"gaus",-1.0,1.0);
      g.SetParameters(h[s]->GetBinContent(h[s]->GetMaximumBin()),0,0.25);
      h[s]->Fit(&g,"RQ0N");
      double mu=g.GetParameter(1); mn=std::min(mn,mu); mx=std::max(mx,mu);
      printf("   %4.0f-%4.0f    %8.0f %10.3f %9.3f\n",lo+s*100.0,lo+(s+1)*100.0,
             h[s]->GetEntries(),mu,2.355*g.GetParameter(2)); }
    printf("   => centroid walks by %.3f MeV across vz 0-1000\n",mx-mn);
  }
}
