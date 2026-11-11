#include "excm.h"
// Where is the (d,d') data in theta_CM, and what Ex does each theta_CM bin cover?
void cmmap(){
  struct C{const char* nm; const char* file; double eb,chi2,keLo,keHi,vzHi,icLo,icHi,kcD,kcP;};
  C cc[2]={
   {"ATTPCROOT","/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/dd/plots_mfi/dd_kin_C15d.root",
    191.25,0.2,0.5,18,1100,1074,1201,3913.0,-1.5},
   {"Spyral","/home/yassid/Data/spyral/C15d_a2091/spyral_kin.root",
    195.0,5.0,0.0,20,400,-1,99999,2535.2,14.5}};
  for(auto&k:cc){
    Pars A; A.ebeam=k.eb; A.chi2=k.chi2; A.thLo=0; A.thHi=180; A.keLo=k.keLo; A.keHi=k.keHi;
    A.vzLo=0; A.vzHi=k.vzHi; A.icLo=k.icLo; A.icHi=k.icHi;
    A.kcDenom=k.kcD; A.kcPivot=k.kcP; A.keOff=0; A.kcOn=true;
    A.ebzSlope=11.0; A.ebzPivot=500; A.ebzOn=true;
    TFile* f=TFile::Open(k.file); if(!f||f->IsZombie()){printf("%s MISSING\n",k.nm);continue;}
    TTree* t=(TTree*)f->Get("dd");
    double ke,th,vz,c2,icd; float icf; int np;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
    t->SetBranchAddress("chi2ndf",&c2);
    bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
    if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
    t->SetBranchAddress("npulse",&np);
    auto* hcm=new TH1D(Form("c%s",k.nm),"",180,0,180); hcm->SetDirectory(nullptr);
    auto* h2=new TH2D(Form("m%s",k.nm),"",36,0,180,120,-2,22); h2->SetDirectory(nullptr);
    long n=0;
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      double icv=icF?icf:icd;
      if(c2>A.chi2||ke<A.keLo||ke>A.keHi||vz<0||vz>A.vzHi||icv<A.icLo||icv>A.icHi) continue;
      double thf=std::fabs(th); if(thf>90) thf=180-thf;
      double thc=thCorr(A,thf,ke);
      double eb=ebAt(A,vz,true);
      double e=exOf(A,eb,thc,ke); if(!std::isfinite(e)) continue;
      double tcm=thcmOf(A,eb,thc,ke);
      hcm->Fill(tcm); h2->Fill(tcm,e); n++; }
    f->Close();
    printf("\n== %s  %ld tracks\n",k.nm,n);
    printf("   theta_CM profile (10 deg bins):\n   ");
    for(int b=0;b<18;b++){ double s=hcm->Integral(hcm->FindBin(b*10+0.5),hcm->FindBin(b*10+9.5));
      if(s>0) printf(" %d-%d:%.1f%%",b*10,b*10+10,100.0*s/n); } printf("\n");
    printf("   mean Ex and rms per 10 deg theta_CM bin:\n");
    for(int b=0;b<18;b++){
      TH1D* py=h2->ProjectionY(Form("p%d",b),h2->GetXaxis()->FindBin(b*10+0.5),h2->GetXaxis()->FindBin(b*10+9.5));
      if(py->GetEntries()<200) continue;
      printf("     %3d-%3d  N %6.0f   <Ex> %6.2f  rms %5.2f   Ex at max %6.2f\n",
             b*10,b*10+10,py->GetEntries(),py->GetMean(),py->GetRMS(),
             py->GetBinCenter(py->GetMaximumBin())); }
  }
}
