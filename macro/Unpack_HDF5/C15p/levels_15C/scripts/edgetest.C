#include "psmodel.h"
#include "exlib.h"
// DECISIVE TEST.  A real state sits at the SAME Ex in every theta slice that can reach it.
// A kinematic edge (the Jacobian piling up at theta_max) MOVES with the slice.
double thmax_of(double Ex,double Kb){
  double Eb=Kb+PS::M15C,pb=std::sqrt(Kb*Kb+2*Kb*PS::M15C);
  TLorentzVector W(0,0,pb,Eb+PS::MD); double s=W.M2(),rs=std::sqrt(s);
  double M4=PS::M15C+Ex,M3=PS::MD; if(rs<=M3+M4) return 0;
  double pcm=PS::omega(s,M3*M3,M4*M4)/(2*rs),E3=std::sqrt(pcm*pcm+M3*M3);
  double bd=pcm/E3,bcm=W.BoostVector().Mag();
  return bd>=bcm?90.0:std::asin(bd/bcm)*180.0/M_PI; }
double ExAtThmax(double th,double Kb){           // invert: which Ex has theta_max = th?
  double lo=0,hi=22.5; for(int i=0;i<60;i++){ double m=0.5*(lo+hi);
    if(thmax_of(m,Kb)>th) lo=m; else hi=m; } return 0.5*(lo+hi); }
void edgetest(){
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
  const int NS=7; double T0=40,DT=4;
  TH1D* h[NS]; for(int i=0;i<NS;i++){h[i]=new TH1D(Form("s%d",i),"",120,0,24);h[i]->SetDirectory(nullptr);}
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
    double icv=icF?icf:icd;
    if(c2>0.2||ke<0.5||ke>18||vz<0||vz>1100||icv<1074||icv>1201) continue;
    double thf=std::fabs(th); if(thf>90) thf=180-thf;
    double thc=thCorr(A,thf,ke);
    double e=exOf(A,ebAt(A,vz,true),thc,ke); if(!std::isfinite(e)) continue;
    int k=(int)((thc-T0)/DT); if(k>=0&&k<NS) h[k]->Fill(e); }
  f->Close();
  printf("For each theta slice: where the counts actually peak, vs the KINEMATIC EDGE Ex(theta_max=theta)\n");
  printf("A real state stays put; an edge tracks the slice.\n\n");
  printf(" slice     N     peak Ex   Ex at kin. edge (slice lo/hi)   counts at Ex=4.65  at Ex=6.87\n");
  for(int i=0;i<NS;i++){
    double tl=T0+i*DT, th2=T0+(i+1)*DT;
    if(h[i]->GetEntries()<500){ printf(" %2.0f-%2.0f  %6.0f   (too few)\n",tl,th2,h[i]->GetEntries()); continue; }
    int bm=0; double mx=0;
    for(int b=1;b<=h[i]->GetNbinsX();b++){ double x=h[i]->GetBinCenter(b);
      if(x<2.0) continue; if(h[i]->GetBinContent(b)>mx){mx=h[i]->GetBinContent(b);bm=b;} }
    printf(" %2.0f-%2.0f  %6.0f    %6.2f        %5.2f  -  %5.2f            %6.0f        %6.0f\n",
      tl,th2,h[i]->GetEntries(),h[i]->GetBinCenter(bm),ExAtThmax(th2,Kb),ExAtThmax(tl,Kb),
      h[i]->GetBinContent(h[i]->FindBin(4.65)), h[i]->GetBinContent(h[i]->FindBin(6.87)));
  }
}
