#include "edgetest.C"
// THE non-circular test.  A real state at Ex* is populated over a RANGE of theta (each theta
// paired with the KE that gives Ex*), ending at theta_max(Ex*).  A Jacobian edge artefact puts
// everything AT theta_max(Ex*).  So: take a narrow Ex window and look at its theta distribution.
void statetest(){
  double Kb=191.25;
  Pars A; A.ebeam=191.25; A.chi2=0.2; A.thLo=0; A.thHi=180; A.keLo=0.5; A.keHi=18;
  A.vzLo=0; A.vzHi=1100; A.icLo=1074; A.icHi=1201;
  A.kcDenom=3913.0; A.kcPivot=-1.5; A.keOff=0; A.kcOn=true;
  A.ebzSlope=11.0; A.ebzPivot=500; A.ebzOn=true;
  TFile* f=TFile::Open("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/dd/plots_mfi/dd_kin_C15d.root");
  TTree* t=(TTree*)f->Get("dd");
  double ke,th,vz,c2,icd; float icf; int np;
  t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
  t->SetBranchAddress("chi2ndf",&c2);
  bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
  if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
  t->SetBranchAddress("npulse",&np);
  struct W{const char* n; double lo,hi;};
  W win[4]={{"state A 4.4-5.0",4.4,5.0},{"off-peak 5.4-6.0",5.4,6.0},
            {"state B 6.6-7.2",6.6,7.2},{"off-peak 8.0-8.6",8.0,8.6}};
  TH1D* h[4]; for(int i=0;i<4;i++){h[i]=new TH1D(Form("w%d",i),"",40,40,80);h[i]->SetDirectory(nullptr);}
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
    double icv=icF?icf:icd;
    if(c2>0.2||ke<0.5||ke>18||vz<0||vz>1100||icv<1074||icv>1201) continue;
    double thf=std::fabs(th); if(thf>90) thf=180-thf;
    double thc=thCorr(A,thf,ke);
    double e=exOf(A,ebAt(A,vz,true),thc,ke); if(!std::isfinite(e)) continue;
    for(int k=0;k<4;k++) if(e>=win[k].lo&&e<win[k].hi) h[k]->Fill(thc); }
  f->Close();
  printf("theta distribution INSIDE a narrow Ex window.\n");
  printf("A real state spreads over theta up to theta_max; an edge artefact sits AT theta_max.\n\n");
  for(int k=0;k<4;k++){
    double tmx=thmax_of(0.5*(win[k].lo+win[k].hi),Kb);
    printf("%-18s  N %5.0f   theta_max(Ex) = %.1f deg\n",win[k].n,h[k]->GetEntries(),tmx);
    printf("   theta:"); for(double x=42;x<=78.01;x+=2) printf("%5.0f",x); printf("\n");
    printf("   N    :"); for(double x=42;x<=78.01;x+=2) printf("%5.0f",h[k]->GetBinContent(h[k]->FindBin(x))); printf("\n");
    // how concentrated is it?  fraction within 2 deg of theta_max
    double near=h[k]->Integral(h[k]->FindBin(tmx-2),h[k]->FindBin(tmx+2));
    printf("   within +-2 deg of theta_max: %.0f / %.0f = %.0f %%\n\n",
           near,h[k]->GetEntries(),100*near/std::max(1.0,h[k]->GetEntries()));
  }
}
