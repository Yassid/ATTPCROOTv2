#include "binacc.h"
#include "binbuild.C"
// The cuts were tightened to suppress a background that lives at LOW Ex.  Fitting only the high-Ex
// region puts that background outside the window, so the cuts can be loosened to buy back
// statistics.  How much is there to gain, and does the high-Ex structure survive?
// Metric: counts in Ex 5-9, and the peak-to-valley of the 6.7 structure (6.4-7.0 over 7.8-8.2),
// which is the strongest thing in that window.
void loosen(){
  for(int ic=0;ic<2;ic++){
    Cfg base=CFG[ic];
    printf("\n================ %s : loosening the selection (Ex 5-9) ================\n",base.nm);
    printf("   %-34s %10s %9s %10s %9s\n","cut","N(5-9)","gain","peak/valley","N(gs)");
    double N0=-1;
    struct V{const char* nm; double chi2,keLo,keHi,icLo,icHi;};
    V vs[8]={{"adopted                      ",base.chi2,base.keLo,base.keHi,base.icLo,base.icHi},
             {"chi2 x2                      ",base.chi2*2,base.keLo,base.keHi,base.icLo,base.icHi},
             {"chi2 x5                      ",base.chi2*5,base.keLo,base.keHi,base.icLo,base.icHi},
             {"chi2 x10                     ",base.chi2*10,base.keLo,base.keHi,base.icLo,base.icHi},
             {"chi2 OFF                     ",1e9,base.keLo,base.keHi,base.icLo,base.icHi},
             {"IC gate OFF                  ",base.chi2,base.keLo,base.keHi,-1,99999},
             {"KE window wide               ",base.chi2,0.0,40.0,base.icLo,base.icHi},
             {"chi2 x5 + IC off + KE wide   ",base.chi2*5,0.0,40.0,-1,99999}};
    for(auto&v:vs){
      Cfg k=base; k.chi2=v.chi2; k.keLo=v.keLo; k.keHi=v.keHi; k.icLo=v.icLo; k.icHi=v.icHi;
      k.vzHi=900;
      Pars A=parsOf(k); A.vzLo=300;
      TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
      double ke,th,vz,c2,icd=0; float icf=0;
      t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
      t->SetBranchAddress("chi2ndf",&c2);
      bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
      if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
      long n59=0,ngs=0,npk=0,nvl=0;
      for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
        double icv=icF?icf:icd;
        if(c2>k.chi2||ke<k.keLo||ke>k.keHi||vz<300||vz>900||icv<k.icLo||icv>k.icHi) continue;
        double thf=std::fabs(th); if(thf>90) thf=180-thf;
        double e=exOf(A,ebAt(A,vz,true),thCorr(A,thf,ke),ke); if(!std::isfinite(e)) continue;
        double tcm=thcmProj(A,ebAt(A,vz,true),thCorr(A,thf,ke),ke);
        if(tcm<10||tcm>60) continue;
        if(e>5&&e<9) n59++;
        if(std::fabs(e)<0.5) ngs++;
        if(e>6.4&&e<7.0) npk++;
        if(e>7.8&&e<8.2) nvl++; }
      f->Close();
      if(N0<0) N0=n59;
      printf("   %-34s %10ld %8.2fx %10.3f %9ld\n",v.nm,n59,N0>0?(double)n59/N0:0,
             nvl>0?(double)npk/nvl*(0.4/0.6):0.0,ngs); }
  }
}
