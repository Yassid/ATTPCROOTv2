#include "excm.h"
#include "binbuild.C"
// How much data is there outside the 10-60 deg I have been using, and where is it dense enough
// to justify FINER bins?  (The statistics are hugely concentrated forward.)
void cmrange(){
  for(auto&k:CFG){
    Pars A=parsOf(k);
    TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
    double ke,th,vz,c2,icd; float icf; int np;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
    t->SetBranchAddress("chi2ndf",&c2);
    bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
    if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
    t->SetBranchAddress("npulse",&np);
    // counts per 5 deg, and how many of them are ABOVE Ex 3 (the resonance region)
    long tot[20]={0}, res[20]={0}; long N=0;
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      double icv=icF?icf:icd;
      if(c2>A.chi2||ke<A.keLo||ke>A.keHi||vz<0||vz>A.vzHi||icv<A.icLo||icv>A.icHi) continue;
      double thf=std::fabs(th); if(thf>90) thf=180-thf;
      double thc=thCorr(A,thf,ke); double eb=ebAt(A,vz,true);
      double e=exOf(A,eb,thc,ke); if(!std::isfinite(e)) continue;
      double tcm=thcmProj(A,eb,thc,ke);
      int b=(int)(tcm/5.0); if(b<0||b>=20) continue;
      tot[b]++; if(e>3.0&&e<9.0) res[b]++; N++; }
    f->Close();
    printf("\n== %s   N = %ld\n",k.nm,N);
    printf("   theta_CM   all counts   in Ex 3-9   (the resonance region)\n");
    for(int b=0;b<16;b++){ if(!tot[b]) continue;
      printf("    %2d-%2d     %8ld     %8ld\n",b*5,b*5+5,tot[b],res[b]); }
  }
}
