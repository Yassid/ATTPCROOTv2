#include "binacc.h"
#include "binbuild.C"
// Is vz the SAME variable in the two codes?  ATTPCROOT's ceiling is 1100 and Spyral's 400, so a
// common VZLO is only meaningful if the coordinate itself is common.  Compare the distributions.
void vzcheck(){
  printf("\n================ vz distribution, gated (d,d') events ================\n");
  printf("   %-10s %9s %8s %8s %8s %8s %8s %8s %8s\n","code","N","min","5%","25%","med","75%","95%","max");
  for(int ic=0;ic<2;ic++){
    Cfg& k=CFG[ic]; Pars A=parsOf(k); A.vzLo=-1e9; A.vzHi=1e9;   // NO vz cut, see the raw range
    TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
    double ke,th,vz,c2; float icf;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
    t->SetBranchAddress("chi2ndf",&c2);t->SetBranchAddress("ic",&icf);
    std::vector<double> v;
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      if(c2>k.chi2||ke<k.keLo||ke>k.keHi||icf<k.icLo||icf>k.icHi) continue;
      v.push_back(vz); }
    f->Close(); std::sort(v.begin(),v.end()); size_t N=v.size();
    auto q=[&](double fr){ size_t j=(size_t)(fr*N); if(j>=N)j=N-1; return v[j]; };
    printf("   %-10s %9zu %8.0f %8.0f %8.0f %8.0f %8.0f %8.0f %8.0f\n",
           k.nm,N,v.front(),q(0.05),q(0.25),q(0.50),q(0.75),q(0.95),v.back());
    // where is the 90-180 spike, and how big, relative to a flat baseline
    long n90=0, nAll=0; for(double x:v){ if(x>0&&x<600) nAll++; if(x>=90&&x<=180) n90++; }
    printf("   %-10s vz in [90,180] = %ld of %ld in [0,600] = %.1f%%  (a flat 90/600 would be 15.0%%)\n",
           "",n90,nAll,100.0*n90/std::max(1L,nAll));
    printf("   %-10s adopted window in binbuild: vzLo=%.0f  vzHi=%.0f\n\n","",0.0,k.vzHi);
  }
}
