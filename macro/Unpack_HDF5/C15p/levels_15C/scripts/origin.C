#include "binacc.h"
#include "binbuild.C"
// Within theta_CM bin 0 the deuteron KE is ~1.4 MeV for BOTH the g.s. (Ex~0) and the excess
// (Ex 2.5-4), so KE cannot separate them.  Compare instead the reconstruction observables:
// vertex z, fit chi2, number of pulses (track length proxy), IC.  An instrumental background
// should look different from the g.s. in at least one of them; a nuclear one should not.
void origin(){
  for(auto&k:CFG){
    Pars A=parsOf(k);
    TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
    double ke,th,vz,c2,icd; float icf; int np;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
    t->SetBranchAddress("chi2ndf",&c2);
    bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
    if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
    t->SetBranchAddress("npulse",&np);
    struct S{const char* n; double lo,hi; long N; double svz,svz2,sc2,sc22,snp,snp2,ske;};
    S w[3]={{"g.s.   Ex -0.5..0.5",-0.5,0.5,0,0,0,0,0,0,0,0},
            {"EXCESS Ex  2.5..4.0", 2.5,4.0,0,0,0,0,0,0,0,0},
            {"states Ex  6.0..7.5", 6.0,7.5,0,0,0,0,0,0,0,0}};
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      double icv=icF?icf:icd;
      if(c2>A.chi2||ke<A.keLo||ke>A.keHi||vz<0||vz>A.vzHi||icv<A.icLo||icv>A.icHi) continue;
      double thf=std::fabs(th); if(thf>90) thf=180-thf;
      double thc=thCorr(A,thf,ke); double eb=ebAt(A,vz,true);
      double e=exOf(A,eb,thc,ke); if(!std::isfinite(e)) continue;
      double tcm=thcmProj(A,eb,thc,ke);
      if(tcm<10||tcm>20) continue;                       // bin 0 only
      for(auto&s:w) if(e>=s.lo&&e<s.hi){
        s.N++; s.svz+=vz; s.svz2+=vz*vz; s.sc2+=c2; s.sc22+=c2*c2;
        s.snp+=np; s.snp2+=np*np; s.ske+=ke; } }
    f->Close();
    printf("\n== %s   theta_CM bin 0 (10-20 deg)\n",k.nm);
    printf("   %-22s %7s %10s %12s %12s %8s\n","window","N","<KE>","<vz> +- rms","<chi2>","<npulse>");
    for(auto&s:w){ if(!s.N) continue;
      double mv=s.svz/s.N, rv=std::sqrt(std::max(0.0,s.svz2/s.N-mv*mv));
      double mc=s.sc2/s.N, mn=s.snp/s.N;
      printf("   %-22s %7ld %10.2f %7.0f +-%4.0f %12.4f %8.1f\n",
             s.n,s.N,s.ske/s.N,mv,rv,mc,mn); }
  }
}
