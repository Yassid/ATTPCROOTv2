#include "exlib.h"
TH1D* fine(const char* file,const Pars&p,const char* nm,long& n,double thMin,double thMax){
  TFile* f=TFile::Open(file); if(!f||f->IsZombie()){n=0;return nullptr;}
  TTree* t=(TTree*)f->Get("dd");
  double ke,th,vz,c2,icd; float icf; int np;
  t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);
  t->SetBranchAddress("vz",&vz);t->SetBranchAddress("chi2ndf",&c2);
  bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
  if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
  t->SetBranchAddress("npulse",&np);
  auto* h=new TH1D(nm,Form("%s;E_{x} [MeV];counts / 25 keV",nm),240,-3,3);
  h->SetDirectory(nullptr); n=0;
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
    double icv=icF?icf:icd, keu=ke+p.keOff;
    if(c2>p.chi2||keu<p.keLo||keu>p.keHi) continue;
    if(vz<p.vzLo||vz>p.vzHi||icv<p.icLo||icv>p.icHi||np<p.npLo||np>p.npHi) continue;
    double thf=std::fabs(th); if(thf>90) thf=180-thf;
    if(thf<thMin||thf>thMax) continue;
    double thc=thCorr(p,thf,keu);
    double e=exOf(p,ebAt(p,vz,true),thc,keu);
    if(std::isfinite(e)){ h->Fill(e); n++; } }
  f->Close(); return h;
}

// (a) float the separation: a real 5/2+ must land at 0.740, not wherever the shape wants it.
// (b) allow the SINGLE peak to be asymmetric (different left/right sigma) -- if the second
//     component was only absorbing a tail, it should lose its significance.
double b2f(double*x,double*p){ double v=p[0]+p[1]*x[0];
  v+=p[2]*std::exp(-0.5*std::pow((x[0]-p[3])/p[4],2));
  v+=p[6]*std::exp(-0.5*std::pow((x[0]-p[3]-p[5])/p[4],2)); return v; }
double asym(double*x,double*p){ // p[2] amp, p[3] mu, p[4] sigL, p[5] sigR
  double s = (x[0]<p[3])? p[4] : p[5];
  return p[0]+p[1]*x[0] + p[2]*std::exp(-0.5*std::pow((x[0]-p[3])/s,2)); }
double asym2(double*x,double*p){ // asymmetric main peak + a 0.740 companion
  double s=(x[0]<p[3])?p[4]:p[5];
  double v=p[0]+p[1]*x[0]+p[2]*std::exp(-0.5*std::pow((x[0]-p[3])/s,2));
  v+=p[6]*std::exp(-0.5*std::pow((x[0]-p[3]-0.740)/p[5],2)); return v; }
void probe(const char* file,const Pars&p,const char* nm,double thMin,double thMax){
  long n; TH1D* h=fine(file,p,Form("P_%s_%.0f",nm,thMin),n,thMin,thMax);
  if(!h||n<500){ printf("%-22s too few\n",nm); return; }
  double lo=-1.8,hi=2.4;
  // (a) separation free
  TF1 ff("ff",b2f,lo,hi,7);
  ff.SetParameters(5,0,h->GetMaximum()*0.9,0,0.3,0.740,h->GetMaximum()*0.1);
  ff.SetParLimits(5,0.15,2.0); ff.SetParLimits(4,0.05,1.2);
  ff.SetParLimits(2,0,1e8); ff.SetParLimits(6,0,1e8);
  h->Fit(&ff,"RQ0N");
  // (b) asymmetric single vs asymmetric + companion
  TF1 fa("fa",asym,lo,hi,6);
  fa.SetParameters(5,0,h->GetMaximum(),0,0.35,0.35);
  fa.SetParLimits(4,0.05,1.5); fa.SetParLimits(5,0.05,1.5);
  h->Fit(&fa,"RQ0N");
  TF1 fa2("fa2",asym2,lo,hi,7);
  for(int i=0;i<6;i++) fa2.SetParameter(i,fa.GetParameter(i));
  fa2.SetParameter(6,h->GetMaximum()*0.05);
  fa2.SetParLimits(4,0.05,1.5); fa2.SetParLimits(5,0.05,1.5); fa2.SetParLimits(6,0,1e8);
  h->Fit(&fa2,"RQ0N");
  printf("%-22s | sep FREE -> %.3f +- %.3f | asym 1pk chi2/ndf %5.2f -> +0.740 companion %5.2f, dchi2 %6.1f | ratio %.3f +- %.3f\n",
    nm, ff.GetParameter(5), ff.GetParError(5),
    fa.GetChisquare()/std::max(1,fa.GetNDF()), fa2.GetChisquare()/std::max(1,fa2.GetNDF()),
    fa.GetChisquare()-fa2.GetChisquare(),
    fa2.GetParameter(2)>0? fa2.GetParameter(6)/fa2.GetParameter(2):0,
    fa2.GetParameter(2)>0? fa2.GetParError(6)/fa2.GetParameter(2):0);
}
void gstest2(){
  Pars A; A.ebeam=191.25; A.chi2=0.2; A.keLo=0.5; A.keHi=18; A.vzLo=0; A.vzHi=1100;
  A.icLo=1074; A.icHi=1201; A.kcDenom=3913.0; A.kcPivot=-1.5; A.kcOn=true;
  A.ebzSlope=11.0; A.ebzPivot=500; A.ebzOn=true;
  Pars S=A; S.ebeam=195; S.chi2=5; S.keLo=0; S.keHi=20; S.vzLo=0; S.vzHi=400;
  S.icLo=-1; S.icHi=99999; S.kcDenom=2535.2; S.kcPivot=14.5;
  const char* fa="/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/dd/plots_mfi/dd_kin_C15d.root";
  const char* fs="/home/yassid/Data/spyral/C15d_a2091/spyral_kin.root";
  probe(fa,A,"ATTPCROOT all",0,180);  probe(fa,A,"ATTPCROOT 65-85",65,85);
  probe(fs,S,"Spyral all",0,180);     probe(fs,S,"Spyral 65-85",65,85);
}
