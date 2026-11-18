#include <Math/SpecFuncMathMore.h>
#include "TLorentzVector.h"
#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TLatex.h"
#include "TLine.h"
#include "TStyle.h"
#include "TMath.h"
#include "TLeaf.h"
#include "TBranch.h"
#include <unordered_map>
#include <vector>
#include <cmath>
#include "binacc.h"
#include "psmodel.h"
#include "binbuild.C"
#include "lineshape.C"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
// CALIBRATE THE STATES: fit ABOVE 5 MeV with the state energies FREE.
// Everything quoted so far came from the (p,p') list transplanted onto our g.s. -- the energies
// were never measured here.  Relative to each code's own g.s. that list gives 5.981 / 6.681 /
// 7.381 / 8.281, while NNDC-like values sit ~0.25 lower.  Let the data say.
//   window 5.0-9.5, low-Ex background excluded by construction
//   continuum  p3*(M) q(M) / Ex^n  x acceptance, n free
//   each state: free energy + free amplitude, width and l FIXED from (p,p')
namespace FE {
   const double SN=1.218, AMU=931.494, HBARC=197.327;
   const double R_FM=1.4*std::pow(15.0,1.0/3.0), MU=AMU*14.0/15.0;
   inline double P_ell(double Ed,int l){ if(Ed<=0) return 0;
     double r=std::sqrt(2*MU*Ed)/HBARC*R_FM;
     double F=r*ROOT::Math::sph_bessel(l,r), G=-r*ROOT::Math::sph_neumann(l,r);
     double d=F*F+G*G; return d>0? r/d : 0; }
   inline double BW(double E,double Er,double Gr,int l){
     double Ed=E-SN, Edr=Er-SN; if(Ed<=0||Edr<=0) return 0;
     double Pr=P_ell(Edr,l); if(Pr<=0) return 0;
     double h=0.5*Gr*P_ell(Ed,l)/Pr;
     return (h/TMath::Pi())/((E-Er)*(E-Er)+h*h); }
   // convolved template, cached on (state, energy rounded to 10 keV)
   struct Key{int s; int e; bool operator==(const Key&o)const{return s==o.s&&e==o.e;} };
   struct KH{ size_t operator()(const Key&k)const{ return (size_t)k.s*100003+k.e; } };
   int NS=4; double GAM[8]={0.300,0.300,0.400,0.500}; int LL[8]={2,2,2,2};
   double SIG=0.40, XLO=4.0, XHI=11.0; int NX=701;
   std::unordered_map<Key,std::vector<double>,KH> CACHE;
   const std::vector<double>& tmpl(int s,double Er){
     Key k{s,(int)std::lround(Er*100)};
     auto it=CACHE.find(k); if(it!=CACHE.end()) return it->second;
     std::vector<double> y(NX,0.0); double dx=(XHI-XLO)/(NX-1);
     const double dE=0.004,E0=SN+1e-6,E1=16.0; int M=(int)((E1-E0)/dE)+1;
     std::vector<double> bw(M);
     for(int j=0;j<M;j++) bw[j]=BW(E0+j*dE,Er,GAM[s],LL[s]);
     double cut=5*SIG, nr=1.0/(std::sqrt(2*TMath::Pi())*SIG);
     for(int i=0;i<NX;i++){ double x=XLO+i*dx,sm=0;
       int jl=std::max(0,(int)((x-cut-E0)/dE)), jh=std::min(M-1,(int)((x+cut-E0)/dE));
       for(int j=jl;j<=jh;j++){ double u=(x-(E0+j*dE))/SIG; sm+=bw[j]*std::exp(-0.5*u*u); }
       y[i]=sm*dE*nr; }
     return CACHE.emplace(k,std::move(y)).first->second; }
   inline double ev(const std::vector<double>&y,double x){
     double t=(x-XLO)/((XHI-XLO)/(NX-1)); int i=(int)t;
     if(i<0||i>=NX-1) return 0; double f=t-i; return y[i]*(1-f)+y[i+1]*f; }
   TH1D* H=nullptr; TGraph* ACC=nullptr; double EB=191.25, GSMU=0.0, exLo=5.0, exHi=9.5;
   double model(double X,const double* p){
     double Xt=X-GSMU;
     double v = std::max(0.0,p[0]) * ((Xt>SN)? PS::dens(Xt,EB)*std::pow(Xt,-std::max(0.0,p[1])) : 0.0)
                                   * (ACC? std::max(0.0,ACC->Eval(Xt)) : 1.0);
     for(int s=0;s<NS;s++) v += std::max(0.0,p[2+2*s]) * ev(tmpl(s,p[3+2*s]),X);
     return v; }
   double fcn(const double* p){ double s=0;
     for(int i=1;i<=H->GetNbinsX();i++){ double x=H->GetBinCenter(i);
       if(x<exLo||x>exHi) continue;
       double n=H->GetBinContent(i), m=model(x,p); if(m<1e-9) m=1e-9;
       s += 2.0*(m-n+(n>0? n*std::log(n/m):0.0)); }
     return s; }
}
void freeE59(double exLo=5.0,double exHi=9.5);
// ---- single-state likelihood SCAN: where does the data actually want structure? ----
void escan(double exLo=5.0,double exHi=9.5);
void escan2(double exLo=3.6,double exHi=9.5);
void freeE59(double exLo,double exHi){
  FE::exLo=exLo; FE::exHi=exHi;
  double seeds[4]={5.98,6.68,7.38,8.28};
  for(int ic=0;ic<2;ic++){
    Cfg k=CFG[ic]; k.vzHi=900; if(ic==0) k.chi2=5.0;     // the loosened ATTPCROOT cut
    Pars A=parsOf(k); A.vzLo=300;
    TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
    double ke,th,vz,c2,icd=0; float icf=0;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
    t->SetBranchAddress("chi2ndf",&c2);
    bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
    if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
    TH1D* h=new TH1D(Form("fe%d",ic),"",240,-2,22); h->SetDirectory(nullptr);
    TH1D* hg=new TH1D(Form("fg%d",ic),"",160,-2,2);  hg->SetDirectory(nullptr);
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      double icv=icF?icf:icd;
      if(c2>k.chi2||ke<k.keLo||ke>k.keHi||vz<300||vz>900||icv<k.icLo||icv>k.icHi) continue;
      double thf=std::fabs(th); if(thf>90) thf=180-thf;
      double eb=ebAt(A,vz,true), thc=thCorr(A,thf,ke);
      double e=exOf(A,eb,thc,ke); if(!std::isfinite(e)) continue;
      double tcm=thcmProj(A,eb,thc,ke); if(tcm<10||tcm>60) continue;
      h->Fill(e); hg->Fill(e); }
    f->Close();
    TF1 gq(Form("gq%d",ic),"gaus",-1.0,1.0);
    gq.SetParameters(hg->GetBinContent(hg->GetMaximumBin()),0,0.25); hg->Fit(&gq,"RQ0N");
    FE::GSMU=gq.GetParameter(1); FE::SIG=gq.GetParameter(2); FE::EB=k.eb;
    FE::CACHE.clear(); FE::H=h;
    Win w{14,86,k.keLo,k.keHi};
    FE::ACC=new TGraph();
    for(int i=0;i<=460;i++){ double X=i*0.05; FE::ACC->SetPoint(i,X,accBin(A,k.eb,X,10,60,w)); }
    int npar=2+2*FE::NS;
    auto* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
    m->SetMaxFunctionCalls(2000000); m->SetTolerance(0.1); m->SetPrintLevel(0); m->SetStrategy(1);
    ROOT::Math::Functor fn(&FE::fcn,npar); m->SetFunction(fn);
    double lvl=h->GetBinContent(h->FindBin(8.5)); if(lvl<1) lvl=1;
    m->SetLimitedVariable(0,"C",lvl*20,lvl*0.1,0,lvl*1e6);
    m->SetLimitedVariable(1,"n",2.0,0.02,0.0,9.0);
    for(int s=0;s<FE::NS;s++){
      m->SetLimitedVariable(2+2*s,Form("A%d",s),lvl*0.5,lvl*0.01,0,lvl*500);
      m->SetLimitedVariable(3+2*s,Form("E%d",s),seeds[s],0.01,seeds[s]-0.45,seeds[s]+0.45); }
    m->Minimize(); m->Hesse();
    const double* p=m->X(); const double* e=m->Errors();
    printf("\n================ %s : states FITTED in Ex %.1f-%.1f ================\n",
           k.nm,exLo,exHi);
    printf("   g.s. mu %+.3f   response sigma %.3f (FWHM %.3f)   continuum n = %.2f\n",
           FE::GSMU,FE::SIG,2.355*FE::SIG,p[1]);
    printf("   %-6s %12s %12s %12s %10s\n","state","Ex (fitted)","above g.s.","(p,p') list","amp");
    for(int s=0;s<FE::NS;s++)
      printf("   %-6d %7.3f+-%.3f %12.3f %12.3f %10.1f\n",
             s+1,p[3+2*s],e[3+2*s],p[3+2*s]-FE::GSMU,seeds[s],p[2+2*s]);
    double c2s=0; int nd=0;
    for(int i=1;i<=h->GetNbinsX();i++){ double x=h->GetBinCenter(i),n=h->GetBinContent(i);
      if(x<exLo||x>exHi||n<5) continue; double mm=FE::model(x,p);
      c2s+=(n-mm)*(n-mm)/n; nd++; }
    printf("   -2lnL %.1f   chi2/ndf %.2f\n",m->MinValue(),c2s/std::max(1,nd-npar));
    delete m;
  }
}

// Fit the CONTINUUM ALONE, then add ONE Breit-Wigner and walk its energy across the window.
// The profile of Delta(-2lnL) is a likelihood "peak finder": it shows where the data put
// structure without ever asking the fit to separate overlapping states, which is what made the
// four-state free fit collapse.
void escan(double exLo,double exHi){
  FE::exLo=exLo; FE::exHi=exHi;
  gStyle->SetOptStat(0);
  TCanvas* c=new TCanvas("c","c",1150,520); c->Divide(2,1);
  int col[2]={kAzure+2,kRed+1};
  for(int ic=0;ic<2;ic++){
    Cfg k=CFG[ic]; k.vzHi=900; if(ic==0) k.chi2=5.0;
    Pars A=parsOf(k); A.vzLo=300;
    TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
    double ke,th,vz,c2,icd=0; float icf=0;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
    t->SetBranchAddress("chi2ndf",&c2);
    bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
    if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
    TH1D* h=new TH1D(Form("sc%d",ic),"",240,-2,22); h->SetDirectory(nullptr);
    TH1D* hg=new TH1D(Form("sg%d",ic),"",160,-2,2);  hg->SetDirectory(nullptr);
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      double icv=icF?icf:icd;
      if(c2>k.chi2||ke<k.keLo||ke>k.keHi||vz<300||vz>900||icv<k.icLo||icv>k.icHi) continue;
      double thf=std::fabs(th); if(thf>90) thf=180-thf;
      double eb=ebAt(A,vz,true), thc=thCorr(A,thf,ke);
      double e=exOf(A,eb,thc,ke); if(!std::isfinite(e)) continue;
      double tcm=thcmProj(A,eb,thc,ke); if(tcm<10||tcm>60) continue;
      h->Fill(e); hg->Fill(e); }
    f->Close();
    TF1 gq(Form("sq%d",ic),"gaus",-1.0,1.0);
    gq.SetParameters(hg->GetBinContent(hg->GetMaximumBin()),0,0.25); hg->Fit(&gq,"RQ0N");
    FE::GSMU=gq.GetParameter(1); FE::SIG=gq.GetParameter(2); FE::EB=k.eb;
    FE::CACHE.clear(); FE::H=h;
    Win w{14,86,k.keLo,k.keHi};
    FE::ACC=new TGraph();
    for(int i=0;i<=460;i++){ double X=i*0.05; FE::ACC->SetPoint(i,X,accBin(A,k.eb,X,10,60,w)); }
    double lvl=h->GetBinContent(h->FindBin(8.5)); if(lvl<1) lvl=1;
    // baseline: continuum only
    FE::NS=0;
    auto run=[&](int npar,double* seed,double* lo,double* hi,const char** nm){
      auto* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
      m->SetMaxFunctionCalls(500000); m->SetTolerance(0.1); m->SetPrintLevel(0);
      ROOT::Math::Functor fn(&FE::fcn,npar); m->SetFunction(fn);
      for(int i=0;i<npar;i++) m->SetLimitedVariable(i,nm[i],seed[i],std::max(1e-6,std::fabs(seed[i])*0.01+1e-4),lo[i],hi[i]);
      m->Minimize(); double v=m->MinValue(); delete m; return v; };
    const char* n0[2]={"C","n"};
    double s0[2]={lvl*20,2.0}, l0[2]={0,0}, h0[2]={lvl*1e6,9};
    double L0=run(2,s0,l0,h0,n0);
    // one BW, energy fixed at each scan point
    FE::NS=1;
    auto* g=new TGraph();
    for(double E=exLo-0.3;E<=exHi+0.3;E+=0.05){
      const char* n1[4]={"C","n","A0","E0"};
      double s1[4]={lvl*20,2.0,lvl*0.5,E}, l1[4]={0,0,0,E-1e-6}, h1[4]={lvl*1e6,9,lvl*500,E+1e-6};
      double L=run(4,s1,l1,h1,n1);
      g->SetPoint(g->GetN(),E-FE::GSMU,L0-L); }
    printf("\n==== %s : single-state scan, Ex %.1f-%.1f  (gsMu %+.3f, FWHM %.3f) ====\n",
           k.nm,exLo,exHi,FE::GSMU,2.355*FE::SIG);
    double best=-1,bE=0;
    for(int i=0;i<g->GetN();i++){ double x,y; g->GetPoint(i,x,y); if(y>best){best=y;bE=x;} }
    printf("   strongest at Ex(above g.s.) = %.3f   Delta(-2lnL) = %.1f\n",bE,best);
    // local maxima above 9
    printf("   local maxima:");
    for(int i=1;i<g->GetN()-1;i++){ double a,b,cc,x,y;
      g->GetPoint(i-1,x,a); g->GetPoint(i,x,y); g->GetPoint(i+1,cc,b);
      if(y>a&&y>b&&y>9) printf("  %.2f (%.0f)",x,y); }
    printf("\n");
    c->cd(ic+1); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetGridy();
    g->SetLineColor(col[ic]); g->SetLineWidth(3); g->SetMarkerStyle(20); g->SetMarkerSize(0.6);
    g->SetMarkerColor(col[ic]);
    g->SetTitle(Form("%s   one BW scanned, continuum free;E_{x} above the g.s. [MeV];#Delta(-2lnL) vs continuum alone",k.nm));
    g->Draw("APL");
    double ref[4]={5.981,6.681,7.381,8.281};
    for(double r:ref){ auto* l=new TLine(r,0,r,best*1.05); l->SetLineStyle(2);
      l->SetLineColor(kGray+2); l->Draw(); }
    auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.030); tx->SetTextColor(kGray+3);
    tx->DrawLatex(0.45,0.86,"dashed: the (p,p') list");
  }
  c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/state_scan_highEx.png");
  printf("\n   wrote state_scan_highEx.png\n");
}

// Two follow-ups in one pass:
//   (a) the SAME single-state scan over a WIDER window, so the strong ~4.4 structure is in view
//       and we can see the whole landscape at once;
//   (b) with the strongest state FIXED at its measured position, scan a SECOND one -- does the
//       data want more than one state, and where?
void escan2(double exLo,double exHi){
  FE::exLo=exLo; FE::exHi=exHi;
  gStyle->SetOptStat(0);
  TCanvas* c=new TCanvas("c2","c2",1150,900); c->Divide(2,2);
  int col[2]={kAzure+2,kRed+1};
  for(int ic=0;ic<2;ic++){
    Cfg k=CFG[ic]; k.vzHi=900; if(ic==0) k.chi2=5.0;
    Pars A=parsOf(k); A.vzLo=300;
    TFile* f=TFile::Open(k.file); TTree* t=(TTree*)f->Get("dd");
    double ke,th,vz,c2,icd=0; float icf=0;
    t->SetBranchAddress("ke",&ke);t->SetBranchAddress("theta",&th);t->SetBranchAddress("vz",&vz);
    t->SetBranchAddress("chi2ndf",&c2);
    bool icF=TString(t->GetBranch("ic")->GetLeaf("ic")->GetTypeName())=="Float_t";
    if(icF) t->SetBranchAddress("ic",&icf); else t->SetBranchAddress("ic",&icd);
    TH1D* h=new TH1D(Form("w%d",ic),"",240,-2,22); h->SetDirectory(nullptr);
    TH1D* hg=new TH1D(Form("wg%d",ic),"",160,-2,2);  hg->SetDirectory(nullptr);
    for(Long64_t i=0;i<t->GetEntries();i++){ t->GetEntry(i);
      double icv=icF?icf:icd;
      if(c2>k.chi2||ke<k.keLo||ke>k.keHi||vz<300||vz>900||icv<k.icLo||icv>k.icHi) continue;
      double thf=std::fabs(th); if(thf>90) thf=180-thf;
      double eb=ebAt(A,vz,true), thc=thCorr(A,thf,ke);
      double e=exOf(A,eb,thc,ke); if(!std::isfinite(e)) continue;
      double tcm=thcmProj(A,eb,thc,ke); if(tcm<10||tcm>60) continue;
      h->Fill(e); hg->Fill(e); }
    f->Close();
    TF1 gq(Form("wq%d",ic),"gaus",-1.0,1.0);
    gq.SetParameters(hg->GetBinContent(hg->GetMaximumBin()),0,0.25); hg->Fit(&gq,"RQ0N");
    FE::GSMU=gq.GetParameter(1); FE::SIG=gq.GetParameter(2); FE::EB=k.eb;
    FE::CACHE.clear(); FE::H=h;
    Win w{14,86,k.keLo,k.keHi};
    FE::ACC=new TGraph();
    for(int i=0;i<=460;i++){ double X=i*0.05; FE::ACC->SetPoint(i,X,accBin(A,k.eb,X,10,60,w)); }
    double lvl=h->GetBinContent(h->FindBin(8.5)); if(lvl<1) lvl=1;
    auto run=[&](int npar,double* seed,double* lo,double* hi,const char** nm){
      auto* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
      m->SetMaxFunctionCalls(500000); m->SetTolerance(0.1); m->SetPrintLevel(0);
      ROOT::Math::Functor fn(&FE::fcn,npar); m->SetFunction(fn);
      for(int i=0;i<npar;i++) m->SetLimitedVariable(i,nm[i],seed[i],
          std::max(1e-6,std::fabs(seed[i])*0.01+1e-4),lo[i],hi[i]);
      m->Minimize(); double v=m->MinValue(); delete m; return v; };
    // --- (a) one state, wide window
    FE::NS=1;
    auto* g1=new TGraph(); double best=-1e9,bE=0;
    for(double E=exLo-0.2;E<=exHi+0.2;E+=0.05){
      const char* n1[4]={"C","n","A0","E0"};
      double s1[4]={lvl*20,2.0,lvl*0.5,E}, l1[4]={0,0,0,E-1e-6}, h1[4]={lvl*1e6,9,lvl*500,E+1e-6};
      double L=run(4,s1,l1,h1,n1);
      g1->SetPoint(g1->GetN(),E-FE::GSMU,-L);
      if(-L>best){best=-L;bE=E-FE::GSMU;} }
    printf("\n==== %s  wide scan %.1f-%.1f  (gsMu %+.3f FWHM %.3f) ====\n",
           k.nm,exLo,exHi,FE::GSMU,2.355*FE::SIG);
    printf("   one state: strongest at %.3f\n",bE);
    printf("   local maxima:");
    for(int i=1;i<g1->GetN()-1;i++){ double x0,y0,x1,y1,x2,y2;
      g1->GetPoint(i-1,x0,y0); g1->GetPoint(i,x1,y1); g1->GetPoint(i+1,x2,y2);
      if(y1>y0&&y1>y2) printf("  %.2f",x1); }
    printf("\n");
    // --- (b) first state pinned, scan a second
    FE::NS=2;
    auto* g2=new TGraph(); double best2=-1e9,bE2=0;
    for(double E=exLo-0.2;E<=exHi+0.2;E+=0.05){
      if(std::fabs(E-(bE+FE::GSMU))<0.25) continue;
      const char* n2[6]={"C","n","A0","E0","A1","E1"};
      double s2[6]={lvl*20,2.0,lvl*0.5,bE+FE::GSMU,lvl*0.2,E};
      double l2[6]={0,0,0,bE+FE::GSMU-1e-6,0,E-1e-6};
      double h2[6]={lvl*1e6,9,lvl*500,bE+FE::GSMU+1e-6,lvl*500,E+1e-6};
      double L=run(6,s2,l2,h2,n2);
      g2->SetPoint(g2->GetN(),E-FE::GSMU,-L);
      if(-L>best2){best2=-L;bE2=E-FE::GSMU;} }
    printf("   second state: best at %.3f   gain over one state = %.1f\n",bE2,best2-best);
    for(int pane=0;pane<2;pane++){
      c->cd(1+ic+2*pane); gPad->SetLeftMargin(0.14); gPad->SetBottomMargin(0.14); gPad->SetGridy();
      TGraph* g=pane? g2 : g1;
      g->SetLineColor(col[ic]); g->SetLineWidth(3);
      g->SetTitle(Form("%s   %s;E_{x} above the g.s. [MeV];-(-2lnL)  [higher = better]",
                       k.nm, pane? Form("SECOND state, first pinned at %.2f",bE)
                                  : "ONE state scanned"));
      g->Draw("AL");
      double ref[6]={4.431,4.931,5.981,6.681,7.381,8.281};
      double ymin=g->GetYaxis()->GetXmin(), ymax=g->GetYaxis()->GetXmax();
      for(double r:ref){ if(r<exLo-0.2||r>exHi+0.2) continue;
        auto* l=new TLine(r,ymin,r,ymax); l->SetLineStyle(2); l->SetLineColor(kGray+2); l->Draw(); }
    }
  }
  c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/state_scan_wide.png");
  printf("\n   wrote state_scan_wide.png\n");
}
