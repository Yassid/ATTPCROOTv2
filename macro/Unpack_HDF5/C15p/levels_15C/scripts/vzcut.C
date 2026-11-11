#include "exlib.h"
// Do the two codes agree on the CONTINUUM SHAPE once each code's own instrumental floor (measured
// below threshold, where no 15C decay can contribute) is removed?  And is any difference just the
// vertex cut -- ATTPCROOT uses vz < 1100 mm, Spyral vz < 400.
void one(const char* file,const Pars&p,const char* nm,double vzHi){
  TFile* f=TFile::Open(file); if(!f||f->IsZombie()){printf("%s MISSING\n",file);return;}
  TTree* t=(TTree*)f->Get("dd");
  double ke,th,vz,c2,icd; float icf; int np;
  t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);
  t->SetBranchAddress("vz",&vz);t->SetBranchAddress("chi2ndf",&c2);
  bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
  if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
  t->SetBranchAddress("npulse",&np);
  auto* h=new TH1D(Form("h%s%.0f",nm,vzHi),"",280,-5,23); h->SetDirectory(nullptr);
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
    double icv=icF?icf:icd, keu=ke+p.keOff;
    if(c2>p.chi2||keu<p.keLo||keu>p.keHi) continue;
    if(vz<p.vzLo||vz>vzHi) continue;
    if(icv<p.icLo||icv>p.icHi||np<p.npLo||np>p.npHi) continue;
    double thf=std::fabs(th); if(thf>90) thf=180-thf;
    double thc=thCorr(p,thf,keu); if(thc<p.thLo||thc>p.thHi) continue;
    double e=exOf(p,ebAt(p,vz,true),thc,keu);
    if(std::isfinite(e)) h->Fill(e); }
  f->Close();
  auto avg=[&](double a,double b){ double s=0; int n=0;
    for(int i=1;i<=h->GetNbinsX();i++){ double x=h->GetBinCenter(i);
      if(x>a&&x<b){s+=h->GetBinContent(i);n++;} } return n?s/n:0.0; };
  double fl=avg(-2.5,-1.5), lo=avg(8,10), hi=avg(18,20);
  printf("  %-12s vz<%-5.0f  floor %6.1f | 8-10 %6.1f | 18-20 %6.1f |"
         " floor-subtracted  %6.1f -> %6.1f   FALL %.2fx\n",
         nm,vzHi,fl,lo,hi,lo-fl,hi-fl,(lo-fl)/std::max(0.1,hi-fl));
}
void vzcut(){
  Pars A; A.ebeam=191.25; A.chi2=0.2; A.thLo=0; A.thHi=180; A.keLo=0.5; A.keHi=18;
  A.vzLo=0; A.vzHi=1100; A.icLo=1074; A.icHi=1201;
  A.kcDenom=3913.0; A.kcPivot=-1.5; A.keOff=0; A.kcOn=true;
  A.ebzSlope=11.0; A.ebzPivot=500; A.ebzOn=true;
  Pars S=A; S.ebeam=195; S.chi2=5; S.keLo=0; S.keHi=20; S.vzLo=0;
  S.icLo=-1; S.icHi=99999; S.kcDenom=2535.2; S.kcPivot=14.5;
  const char* FA="/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/dd/plots_mfi/dd_kin_C15d.root";
  const char* FS="/home/yassid/Data/spyral/C15d_a2091/spyral_kin.root";
  printf("continuum level per 100 keV bin; 'floor' is Ex -2.5..-1.5, below threshold\n\n");
  one(FA,A,"ATTPCROOT",1100); one(FA,A,"ATTPCROOT",400);
  one(FS,S,"Spyral",400);     one(FS,S,"Spyral",1100);
}
