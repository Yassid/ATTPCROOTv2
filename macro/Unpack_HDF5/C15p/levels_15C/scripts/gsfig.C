#include "exlib.h"
TH1D* fineH(const char* file,const Pars&p,const char* nm,long& n,double thMin,double thMax){
  TFile* f=TFile::Open(file); if(!f||f->IsZombie()){n=0;return nullptr;}
  TTree* t=(TTree*)f->Get("dd");
  double ke,th,vz,c2,icd; float icf; int np;
  t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);
  t->SetBranchAddress("vz",&vz);t->SetBranchAddress("chi2ndf",&c2);
  bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
  if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
  t->SetBranchAddress("npulse",&np);
  auto* h=new TH1D(nm,"",240,-3,3); h->SetDirectory(nullptr); n=0;
  for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
    double icv=icF?icf:icd, keu=ke+p.keOff;
    if(c2>p.chi2||keu<p.keLo||keu>p.keHi) continue;
    if(vz<p.vzLo||vz>p.vzHi||icv<p.icLo||icv>p.icHi||np<p.npLo||np>p.npHi) continue;
    double thf=std::fabs(th); if(thf>90) thf=180-thf;
    if(thf<thMin||thf>thMax) continue;
    double e=exOf(p,ebAt(p,vz,true),thCorr(p,thf,keu),keu);
    if(std::isfinite(e)){ h->Fill(e); n++; } }
  f->Close(); return h; }
double asym2(double*x,double*p){
  double s=(x[0]<p[3])?p[4]:p[5];
  double v=p[0]+p[1]*x[0]+p[2]*std::exp(-0.5*std::pow((x[0]-p[3])/s,2));
  v+=p[6]*std::exp(-0.5*std::pow((x[0]-p[3]-0.740)/p[5],2)); return v; }
void gsfig(){
  gStyle->SetOptStat(0);
  Pars A; A.ebeam=191.25; A.chi2=0.2; A.keLo=0.5; A.keHi=18; A.vzLo=0; A.vzHi=1100;
  A.icLo=1074; A.icHi=1201; A.kcDenom=3913.0; A.kcPivot=-1.5; A.kcOn=true;
  A.ebzSlope=11.0; A.ebzPivot=500; A.ebzOn=true;
  Pars S=A; S.ebeam=195; S.chi2=5; S.keLo=0; S.keHi=20; S.vzLo=0; S.vzHi=400;
  S.icLo=-1; S.icHi=99999; S.kcDenom=2535.2; S.kcPivot=14.5;
  const char* ff[2]={"/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15d/dd/plots_mfi/dd_kin_C15d.root",
                     "/home/yassid/Data/spyral/C15d_a2091/spyral_kin.root"};
  Pars pp[2]={A,S}; const char* nm[2]={"ATTPCROOT","Spyral"};
  auto* c=new TCanvas("c","",1400,800); c->Divide(1,2);
  for(int k=0;k<2;k++){
    long n; TH1D* h=fineH(ff[k],pp[k],Form("G%d",k),n,0,180);
    TF1 f("f",asym2,-1.8,2.4,7);
    f.SetParameters(5,0,h->GetMaximum(),0,0.35,0.35,h->GetMaximum()*0.05);
    f.SetParLimits(4,0.05,1.5); f.SetParLimits(5,0.05,1.5); f.SetParLimits(6,0,1e8);
    h->Fit(&f,"RQ0N");
    double bw=h->GetBinWidth(1);
    double Ygs=f.GetParameter(2)*0.5*(f.GetParameter(4)+f.GetParameter(5))*std::sqrt(2*M_PI)/bw;
    double Yex=f.GetParameter(6)*f.GetParameter(5)*std::sqrt(2*M_PI)/bw;
    double eYex=f.GetParError(6)*f.GetParameter(5)*std::sqrt(2*M_PI)/bw;
    printf("%-10s n=%6ld | g.s. yield %7.0f | 5/2+ yield %6.0f +- %.0f | ratio %.4f +- %.4f | FWHM(L,R) %.3f %.3f\n",
      nm[k],n,Ygs,Yex,eYex,Yex/Ygs,eYex/Ygs,2.355*f.GetParameter(4),2.355*f.GetParameter(5));
    c->cd(k+1); gPad->SetLogy();
    h->SetTitle(Form("%s  15C(d,d')  g.s. and the 5/2^{+} at 0.740;E_{x} [MeV];counts / 25 keV",nm[k]));
    h->GetXaxis()->SetRangeUser(-2,2.5); h->SetLineColor(kBlack); h->SetLineWidth(2); h->Draw("hist");
    f.SetLineColor(kRed); f.SetLineWidth(2); f.DrawCopy("same");
    auto* g2=new TF1(Form("g2_%d",k),"[0]*exp(-0.5*((x-0.740)/[1])^2)",-1.8,2.4);
    g2->SetParameters(f.GetParameter(6),f.GetParameter(5));
    g2->SetLineColor(kBlue); g2->SetLineStyle(2); g2->SetLineWidth(2); g2->DrawCopy("same");
    auto* l=new TLine(0.740,0.5,0.740,h->GetMaximum()*0.3);
    l->SetLineColor(kGreen+2); l->SetLineStyle(3); l->SetLineWidth(2); l->Draw();
    auto* t=new TLatex(); t->SetTextSize(0.05); t->SetTextColor(kBlue);
    t->DrawLatex(1.2,h->GetMaximum()*0.08,Form("5/2^{+}/g.s. = %.3f",Yex/Ygs));
  }
  c->SaveAs("/home/yassid/C15d_gs_and_5half.png");
  printf("wrote /home/yassid/C15d_gs_and_5half.png\n");
}
