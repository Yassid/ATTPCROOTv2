#include "psmodel.h"
#include "exlib.h"
// Where do 3-body d+14C+n deuterons actually GO, and does the detector see them?
// Compare the simulated (theta, KE) distribution against the region the real data occupies.
void acccheck(){
  Pars A; A.ebeam=191.25; A.keLo=0.5; A.keHi=18; A.thLo=0; A.thHi=180;
  A.kcDenom=3913.0; A.kcPivot=-1.5; A.kcOn=true;

  // --- real data coverage in (theta_corr, KE) ---
  TFile* f=TFile::Open("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/dd/plots_mfi/dd_kin_C15d.root");
  TTree* t=(TTree*)f->Get("dd");
  double ke,th,vz,c2,icd; float icf; int np;
  t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
  t->SetBranchAddress("chi2ndf",&c2);
  bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
  if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
  t->SetBranchAddress("npulse",&np);
  auto* hD=new TH2D("hD","data;theta;KE",90,0,90,80,0,20); hD->SetDirectory(nullptr);
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
    double icv=icF?icf:icd;
    if(c2>0.2||ke<0.5||ke>18||vz<0||vz>1100||icv<1074||icv>1201) continue;
    double thf=std::fabs(th); if(thf>90) thf=180-thf;
    hD->Fill(thCorr(A,thf,ke),ke); }
  f->Close();
  printf("data: %.0f tracks in the (theta,KE) map\n",hD->GetEntries());
  printf("  theta range holding 99%% of them: ");
  { TH1D* px=hD->ProjectionX(); double tot=px->Integral(),c=0; double lo=-1,hi=-1;
    for(int i=1;i<=px->GetNbinsX();i++){ c+=px->GetBinContent(i);
      if(lo<0&&c>0.005*tot) lo=px->GetBinCenter(i);
      if(hi<0&&c>0.995*tot) hi=px->GetBinCenter(i); }
    printf("%.0f - %.0f deg\n",lo,hi); }

  // --- simulated 3-body: where do the deuterons go? ---
  double Kb=191.25;
  TLorentzVector beam(0,0,std::sqrt(Kb*Kb+2*Kb*PS::M15C),Kb+PS::M15C);
  TLorentzVector W=beam+TLorentzVector(0,0,0,PS::MD);
  double m3[3]={PS::MD,PS::M14C,PS::MN}; TGenPhaseSpace g; g.SetDecay(W,3,m3);
  auto* hS=new TH2D("hS","sim;theta;KE",90,0,90,80,0,20); hS->SetDirectory(nullptr);
  double wAll=0,wKE=0,wTh=0,wBoth=0;
  for(int i=0;i<2000000;i++){ double w=g.Generate();
    TLorentzVector* d=g.GetDecay(0);
    double k=d->E()-PS::MD, th0=d->Theta()*180.0/M_PI;
    wAll+=w;
    bool okKE=(k>=0.5&&k<=18), okTh=(th0>=40&&th0<=90);
    if(okKE) wKE+=w; if(okTh) wTh+=w; if(okKE&&okTh){ wBoth+=w; hS->Fill(th0,k,w); } }
  printf("\nsimulated 3-body deuterons (phase-space weighted):\n");
  printf("  pass KE 0.5-18            : %5.1f %%\n",100*wKE/wAll);
  printf("  pass theta 40-90 deg      : %5.1f %%\n",100*wTh/wAll);
  printf("  pass BOTH                 : %5.1f %%   <- what the detector could ever see\n",100*wBoth/wAll);
  printf("\n  simulated theta profile (phase-space weight, all KE):\n   ");
  { TGenPhaseSpace g2; g2.SetDecay(W,3,m3);
    double prof[9]={0}; double tot=0;
    for(int i=0;i<1000000;i++){ double w=g2.Generate(); TLorentzVector* d=g2.GetDecay(0);
      double th0=d->Theta()*180.0/M_PI; int b=(int)(th0/10); if(b>8)b=8; prof[b]+=w; tot+=w; }
    for(int b=0;b<9;b++) printf(" %d-%d:%4.1f%%",b*10,b*10+10,100*prof[b]/tot); printf("\n"); }
  printf("\n  data theta profile:\n   ");
  { TH1D* px=hD->ProjectionX(); double tot=px->Integral();
    for(int b=0;b<9;b++){ double s=px->Integral(px->FindBin(b*10+0.5),px->FindBin(b*10+9.5));
      printf(" %d-%d:%4.1f%%",b*10,b*10+10,100*s/tot); } printf("\n"); }
}
