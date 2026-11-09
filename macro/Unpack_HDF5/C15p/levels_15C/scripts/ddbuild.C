#include "exlib.h"
// Build the (d,d') Ex spectrum with the user's SAVED explorer parameters, reproducing the page's
// own transforms (kc theta-correction and Ebeam(vz)) so the offline numbers match the screen.
TH1D* build(const char* file,const Pars&p,const char* nm,long& nOut,bool thFold=true){
  TFile* f=TFile::Open(file); if(!f||f->IsZombie()){printf("%s MISSING\n",file);nOut=0;return nullptr;}
  TTree* t=(TTree*)f->Get("dd");
  double ke,th,vz,c2,icd; float icf; int np;
  t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);
  t->SetBranchAddress("vz",&vz);t->SetBranchAddress("chi2ndf",&c2);
  bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
  if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
  t->SetBranchAddress("npulse",&np);
  auto* h=new TH1D(nm,Form("%s;E_{x} [MeV];counts / 100 keV",nm),280,-5,23);
  h->SetDirectory(nullptr);
  long n=0;
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
    double icv=icF?icf:icd;
    double keu=ke+p.keOff;
    if(c2>p.chi2) continue;
    if(keu<p.keLo||keu>p.keHi) continue;
    if(vz<p.vzLo||vz>p.vzHi) continue;
    if(icv<p.icLo||icv>p.icHi) continue;
    if(np<p.npLo||np>p.npHi) continue;
    double thf=std::fabs(th); if(thFold&&thf>90) thf=180-thf;
    double thc=thCorr(p,thf,keu);
    if(thc<p.thLo||thc>p.thHi) continue;
    double eb=ebAt(p,vz,true);
    double e=exOf(p,eb,thc,keu);
    if(std::isfinite(e)){ h->Fill(e); n++; } }
  f->Close(); nOut=n; return h;
}
void ddbuild(){
  gStyle->SetOptStat(0);
  // ---- ATTPCROOT, from explorer_params_..._atttpcroot.json
  Pars A; A.ebeam=191.25; A.chi2=0.2; A.thLo=0; A.thHi=180; A.keLo=0.5; A.keHi=18;
  A.vzLo=0; A.vzHi=1100; A.icLo=1074; A.icHi=1201;
  A.kcDenom=3913.0; A.kcPivot=-1.5; A.keOff=0; A.kcOn=true;
  A.ebzSlope=11.0; A.ebzPivot=500; A.ebzOn=true;
  // ---- Spyral, from explorer_params_..._spyral.json
  Pars S=A; S.ebeam=195; S.chi2=5; S.keLo=0; S.keHi=20; S.vzLo=0; S.vzHi=400;
  S.icLo=-1; S.icHi=99999; S.kcDenom=2535.2; S.kcPivot=14.5;
  long na,ns;
  TH1D* ha=build("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/dd/plots_mfi/dd_kin_C15d.root",A,"ATTPCROOT",na);
  TH1D* hs=build("/home/yassid/Data/spyral/C15d_a2091/spyral_kin.root",S,"Spyral",ns);
  TFile o("/home/yassid/Data/C15d_statefit/dd_ex.root","RECREATE");
  for(auto* h : {ha,hs}){ if(!h) continue;
    TF1 g("g","gaus(0)+pol1(3)",-2.0,2.0);
    g.SetParameters(h->GetMaximum(),0,0.6,10,0); g.SetParLimits(2,0.08,3.0);
    h->Fit(&g,"RQ0N");
    printf("%-10s %7ld tracks | g.s. mu %+.3f  FWHM %.3f\n",
           h->GetName(), (h==ha?na:ns), g.GetParameter(1), 2.355*fabs(g.GetParameter(2)));
    h->Write(h->GetName()); }
  o.Close();
  printf("wrote dd_ex.root\n");
}
