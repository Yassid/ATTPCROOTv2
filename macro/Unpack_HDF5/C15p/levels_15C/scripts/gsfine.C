#include "exlib.h"
// Fine-binned elastic region, to test whether the 5/2+ at 0.740 is required alongside the g.s.
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
// models
double b1(double*x,double*p){ // single peak + linear
  return p[0]+p[1]*x[0] + p[2]*std::exp(-0.5*std::pow((x[0]-p[3])/p[4],2)); }
double b2(double*x,double*p){ // g.s. + 5/2+ at a fixed separation p[5], common sigma
  double v=p[0]+p[1]*x[0];
  v += p[2]*std::exp(-0.5*std::pow((x[0]-p[3])/p[4],2));
  v += p[6]*std::exp(-0.5*std::pow((x[0]-p[3]-p[5])/p[4],2));
  return v; }
void test(TH1D* h,const char* nm,long n){
  if(!h||n<500){ printf("%-28s too few\n",nm); return; }
  double lo=-1.8,hi=2.4;
  TF1 f1("f1",b1,lo,hi,5);
  f1.SetParameters(5,0,h->GetMaximum(),0,0.35); f1.SetParLimits(4,0.05,1.2);
  h->Fit(&f1,"RQ0N");
  TF1 f2("f2",b2,lo,hi,7);
  f2.SetParameters(f1.GetParameter(0),f1.GetParameter(1),f1.GetParameter(2)*0.9,
                   f1.GetParameter(3),f1.GetParameter(4)*0.8,0.740,f1.GetParameter(2)*0.15);
  f2.FixParameter(5,0.740);
  f2.SetParLimits(4,0.05,1.2); f2.SetParLimits(2,0,1e8); f2.SetParLimits(6,0,1e8);
  h->Fit(&f2,"RQ0N");
  double c1=f1.GetChisquare(), c2=f2.GetChisquare();
  int d1=f1.GetNDF(), d2=f2.GetNDF();
  double A=f2.GetParameter(2), B=f2.GetParameter(6), eB=f2.GetParError(6), s=fabs(f2.GetParameter(4));
  printf("%-28s n=%6ld | 1 peak chi2/ndf %6.2f | 2 peaks %6.2f | dchi2 %7.1f (1 dof) | FWHM %.3f | 5/2+ / g.s. = %.3f +- %.3f\n",
         nm, n, c1/std::max(1,d1), c2/std::max(1,d2), c1-c2, 2.355*s, A>0?B/A:0, A>0?eB/A:0);
}
void gsfine(){
  Pars A; A.ebeam=191.25; A.chi2=0.2; A.keLo=0.5; A.keHi=18; A.vzLo=0; A.vzHi=1100;
  A.icLo=1074; A.icHi=1201; A.kcDenom=3913.0; A.kcPivot=-1.5; A.kcOn=true;
  A.ebzSlope=11.0; A.ebzPivot=500; A.ebzOn=true;
  Pars S=A; S.ebeam=195; S.chi2=5; S.keLo=0; S.keHi=20; S.vzLo=0; S.vzHi=400;
  S.icLo=-1; S.icHi=99999; S.kcDenom=2535.2; S.kcPivot=14.5;
  const char* fa="/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/dd/plots_mfi/dd_kin_C15d.root";
  const char* fs="/home/yassid/Data/spyral/C15d_a2091/spyral_kin.root";
  long n;
  struct W{double lo,hi;const char*l;};
  W ws[4]={{0,180,"all theta"},{55,75,"theta 55-75"},{65,85,"theta 65-85"},{70,88,"theta 70-88"}};
  for(auto&w:ws){ TH1D* h=fine(fa,A,Form("A_%s",w.l),n,w.lo,w.hi); test(h,Form("ATTPCROOT %s",w.l),n); }
  printf("\n");
  for(auto&w:ws){ TH1D* h=fine(fs,S,Form("S_%s",w.l),n,w.lo,w.hi); test(h,Form("Spyral %s",w.l),n); }
}
