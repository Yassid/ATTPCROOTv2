#include "edgetest.C"
// Counts in the KINEMATICALLY FORBIDDEN region (Ex above theta_max for that track's own theta)
// cannot be 15C(d,d').  They measure the background directly, with no model.
void forbid(){
  double Kb=191.25;
  struct Cfg{const char* nm; const char* file; double eb,chi2,keLo,keHi,vzHi,icLo,icHi,kcD,kcP;};
  Cfg C[2]={
   {"ATTPCROOT","/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/dd/plots_mfi/dd_kin_C15d.root",
    191.25,0.2,0.5,18,1100,1074,1201,3913.0,-1.5},
   {"Spyral","/home/yassid/Data/spyral/C15d_a2091/spyral_kin.root",
    195.0,5.0,0.0,20,400,-1,99999,2535.2,14.5}};
  for(auto&c:C){
    Pars A; A.ebeam=c.eb; A.chi2=c.chi2; A.thLo=0; A.thHi=180; A.keLo=c.keLo; A.keHi=c.keHi;
    A.vzLo=0; A.vzHi=c.vzHi; A.icLo=c.icLo; A.icHi=c.icHi;
    A.kcDenom=c.kcD; A.kcPivot=c.kcP; A.keOff=0; A.kcOn=true;
    A.ebzSlope=11.0; A.ebzPivot=500; A.ebzOn=true;
    TFile* f=TFile::Open(c.file); if(!f||f->IsZombie()){printf("%s MISSING\n",c.nm);continue;}
    TTree* t=(TTree*)f->Get("dd");
    double ke,th,vz,c2,icd; float icf; int np;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
    t->SetBranchAddress("chi2ndf",&c2);
    bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
    if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
    t->SetBranchAddress("npulse",&np);
    long tot=0, forb=0, forbHi=0;
    // Ex spectrum of the ALLOWED and FORBIDDEN samples
    auto* hA=new TH1D(Form("A%s",c.nm),"",280,-5,23); hA->SetDirectory(nullptr);
    auto* hF=new TH1D(Form("F%s",c.nm),"",280,-5,23); hF->SetDirectory(nullptr);
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      double icv=icF?icf:icd;
      if(c2>A.chi2||ke<A.keLo||ke>A.keHi||vz<0||vz>A.vzHi||icv<A.icLo||icv>A.icHi) continue;
      double thf=std::fabs(th); if(thf>90) thf=180-thf;
      double thc=thCorr(A,thf,ke);
      double e=exOf(A,ebAt(A,vz,true),thc,ke); if(!std::isfinite(e)) continue;
      tot++;
      double tmx=thmax_of(e,c.eb);        // largest lab angle this Ex can produce
      if(thc>tmx+1.0){ forb++; hF->Fill(e); if(e>3) forbHi++; }
      else hA->Fill(e);
    }
    f->Close();
    printf("\n== %s  (theta corrected; 'forbidden' = theta exceeds theta_max(Ex) by >1 deg)\n",c.nm);
    printf("   total %ld | forbidden %ld (%.1f %%) | forbidden with Ex>3 %ld (%.1f %%)\n",
           tot,forb,100.0*forb/tot,forbHi,100.0*forbHi/tot);
    printf("   Ex:        "); for(double x=0;x<=18.01;x+=2) printf("%7.0f",x); printf("\n");
    printf("   allowed:   "); for(double x=0;x<=18.01;x+=2) printf("%7.0f",hA->GetBinContent(hA->FindBin(x))); printf("\n");
    printf("   FORBIDDEN: "); for(double x=0;x<=18.01;x+=2) printf("%7.0f",hF->GetBinContent(hF->FindBin(x))); printf("\n");
  }
}
