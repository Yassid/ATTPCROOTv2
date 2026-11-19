#include "resfitD.C"
// ONE ANGULAR BIN (10-15), fitted over the FULL range: g.s. through 10 MeV, with NO gap.
// Fitting only 3.5-9 meant the background was fitted on the falling tail alone, which gave a
// steep power law that is nonsense below 3.5 -- it diverged at threshold.  Here the background
// must span everything, so it has to turn ON at Sn and describe the rise to ~4 as well as the
// fall beyond.
//   free: g.s. amplitude, 5/2+ amplitude, background normalisation + shape, 4 group amplitudes
//   fixed: the response (from the g.s. peak), the state energies, the pair ratios
namespace B1 {
  int BIN=0, FORM=0, NSH=1;
  TH1D* H=nullptr; TGraph* ACC=nullptr; double EB=191.25;
  double exLo=-1.5, exHi=10.0;
  const double SN=PS::SN;
  // candidate background shapes, all zero below Sn
  double shp(double Ex,const double* s){
    if(Ex<=SN) return 0.0;
    double d=PS::dens(Ex,EB); if(d<=0) return 0.0;
    double u=Ex-SN;
    switch(FORM){
      case 0: return d;                                   // pure phase space, 0 shape par
      case 1: return d*std::pow(Ex,-std::max(0.0,s[0]));  // x 1/Ex^n
      case 2: return d*std::exp(-std::max(0.0,s[0])*u);   // x exp(-k(Ex-Sn))
      case 3: return d*std::pow(u,std::max(0.0,s[0]))*std::exp(-std::max(0.0,s[1])*u);
      case 4: { double w=std::max(0.05,s[0]);             // smooth turn-on x exponential fall
                return (1.0-std::exp(-u/w))*std::exp(-std::max(0.0,s[1])*Ex); }
      case 5: { double w=std::max(0.05,s[0]);             // turn-on x phase space x fall
                return d*(1.0-std::exp(-u/w))*std::exp(-std::max(0.0,s[1])*u); }
    }
    return d;
  }
  // ---- g.s. and 5/2+ as PURE GAUSSIANS.
  // They were a bifurcated Gaussian with exponential tails (shapeGT), whose widths and tail
  // onsets came from the elastic peak fit and were then EXTRAPOLATED under the resonances --
  // never checked.  Here the centroid and a single COMMON width are fitted instead, and they are
  // the last two parameters of the vector.
  bool GAUSS=true;        // kept only so old call sites compile; the peaks are ALWAYS Gaussian
  int  NPAR=0;
  double PKHI=2.5;        // |Ex| below this is the "peak region", above it the physics region
  bool PEAKFREE=true;     // false: the peak is outside the fit, so no g.s./5/2+ parameters
  double MUFIX=0.0;       // the g.s. centroid, measured separately, used when PEAKFREE is false
  double model(double X,const double* p){
    RK::Bin& d=RK::B[BIN];
    double v=0, mu=MUFIX;
    if(PEAKFREE){
      mu = p[NPAR-2];
      double sig = std::max(0.05,p[NPAR-1]);
      double t0=(X-mu)/sig, t1=(X-mu-0.740)/sig;
      v = std::max(0.0,p[0])*std::exp(-0.5*t0*t0)
        + std::max(0.0,p[1])*std::exp(-0.5*t1*t1); }
    double Xt=X-mu;
    const int o = PEAKFREE ? 2 : 0;
    v += std::max(0.0,p[o])*shp(Xt,p+o+1)*(ACC?std::max(0.0,ACC->Eval(Xt)):1.0);
    const double* a=p+o+1+NSH;
    for(int r=0;r<RK::NR;r++) if(RK::GRP[r]>=0&&d.in[r]){
      double amp=a[RK::GRP[r]]*a[RK::GRP[r]];          // s^2 : >= 0 and smooth through 0
      v += amp*((RK::RIDX[r]>=0)?RK::FIXRAT[RK::RIDX[r]]:1.0)*d.T[r].Eval(X); }
    return v;
  }
  double fcn(const double* p){
    RK::Bin& d=RK::B[BIN]; double s=0;
    for(int i=1;i<=d.h->GetNbinsX();i++){
      double x=d.h->GetBinCenter(i);
      if(x<exLo||x>exHi) continue;
      double n=d.h->GetBinContent(i), m=model(x,p); if(m<1e-9) m=1e-9;
      s += 2.0*(m-n+(n>0? n*std::log(n/m):0.0)); }
    return s; }
}
void one1015(int ic=0,int drawForm=3);
void one1015(int ic,int drawForm){
  setCM(10,5,60);
  RK::exLo=3.5; RK::exHi=9.0; RK::ONETERM=true; RK::BKG=RK::HALON;
  RK::TIE12=false; RK::PAIRS=true; RK::DELTA=0.0;
  for(int i=0;i<RK::NR;i++) RK::USE[i]=(i<6);
  Cfg b0=CFG[0], b1=CFG[1];
  CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; if(ic==0) CFG[ic].chi2=5.0;
  VZLO=300; BB_QUIET=true; binbuild();
  RK::SEED.clear(); rkRun(CFG[ic],NCM,false);
  for(int r=0;r<RK::NRAT;r++) RK::FIXRAT[r]=RK::LASTP[1+NCM*RK::nper()+r];
  B1::BIN=0; B1::EB=CFG[ic].eb;
  RK::Bin& d=RK::B[0];
  Pars A=parsOf(CFG[ic]); A.vzLo=300; Win w{14,86,CFG[ic].keLo,CFG[ic].keHi};
  B1::ACC=new TGraph();
  for(int i=0;i<=460;i++){ double X=i*0.05; B1::ACC->SetPoint(i,X,accBin(A,CFG[ic].eb,X,d.cmLo,d.cmHi,w)); }
  struct F{int id; int nsh; const char* nm;};
  F fs[6]={{0,0,"p3* q                      pure phase space  "},
           {1,1,"p3* q / Ex^n                                 "},
           {2,1,"p3* q exp(-k(Ex-Sn))                         "},
           {3,2,"p3* q (Ex-Sn)^a exp(-k(Ex-Sn))               "},
           {4,2,"(1-exp(-(Ex-Sn)/w)) exp(-k Ex)               "},
           {5,2,"p3* q (1-exp(-(Ex-Sn)/w)) exp(-k(Ex-Sn))     "}};
  printf("\n================ %s  theta_CM %.0f-%.0f : background over the FULL range %.1f to %.1f ================\n",
         CFG[ic].nm,d.cmLo,d.cmHi,B1::exLo,B1::exHi);
  printf("   g.s. mu %+.3f  FWHM %.2f/%.2f   pair ratios %.2f/%.2f\n",
         d.gsMu,2.355*d.sgL,2.355*d.sgR,RK::FIXRAT[0],RK::FIXRAT[1]);
  printf("   %-46s %4s %10s %9s %s\n","background form","npar","-2lnL","chi2/ndf","shape");
  for(auto&f:fs){
    B1::FORM=f.id; B1::NSH=f.nsh;
    int npar=3+f.nsh+4;
    auto* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
    m->SetMaxFunctionCalls(500000); m->SetTolerance(0.01); m->SetPrintLevel(0); m->SetStrategy(2);
    ROOT::Math::Functor fn(&B1::fcn,npar); m->SetFunction(fn);
    double pk=d.h->GetBinContent(d.h->GetMaximumBin());
    m->SetLimitedVariable(0,"Ags",pk,pk*0.01,0,pk*10);
    m->SetLimitedVariable(1,"A52",0.1*pk,pk*0.001,0,pk*5);
    m->SetLimitedVariable(2,"B",50,1,0,1e7);
    for(int s=0;s<f.nsh;s++) m->SetLimitedVariable(3+s,Form("s%d",s),
        (f.id==3&&s==0)?0.5:((f.id>=4&&s==0)?1.0:0.5),0.02,0.0,(f.id>=4&&s==0)?8.0:9.0);
    for(int g=0;g<4;g++) m->SetLimitedVariable(3+f.nsh+g,Form("A%d",g),20,0.5,0,1e5);
    m->Minimize();
    const double* p=m->X();
    double s2=0; int nd=0;
    for(int i=1;i<=d.h->GetNbinsX();i++){ double x=d.h->GetBinCenter(i),n=d.h->GetBinContent(i);
      if(x<B1::exLo||x>B1::exHi) continue;
      double mm=B1::model(x,p); if(mm<1e-9) mm=1e-9; s2+=(n-mm)*(n-mm)/std::max(mm,1.0); nd++; }
    TString sh;
    for(int s=0;s<f.nsh;s++) sh+=Form("%6.3f ",p[3+s]);
    printf("   %-46s %4d %10.1f %9.2f  %s\n",f.nm,npar,m->MinValue(),s2/std::max(1.0,nd-(double)npar),sh.Data());
    if(f.id==drawForm){
      std::vector<double> bp(p,p+npar);
      gStyle->SetOptStat(0);
      TCanvas* c=new TCanvas("c1","c1",900,640);
      gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.13); gPad->SetLogy();
      auto* hd=(TH1D*)d.h->Clone("hd1015"); hd->SetDirectory(nullptr);
      hd->GetXaxis()->SetRangeUser(B1::exLo,B1::exHi);
      double mx=0; for(int i=1;i<=hd->GetNbinsX();i++){ double x=hd->GetBinCenter(i);
        if(x>=B1::exLo&&x<=B1::exHi) mx=std::max(mx,hd->GetBinContent(i)); }
      hd->SetMaximum(4*std::max(1.0,mx)); hd->SetMinimum(0.4);
      hd->SetMarkerStyle(20); hd->SetMarkerSize(0.6); hd->SetLineColor(kBlack);
      hd->SetTitle(Form("%s   #theta_{CM} %.0f-%.0f#circ   fitted g.s. to 10 MeV, NO gap;E_{x} [MeV];counts / 100 keV",
                        CFG[ic].nm,d.cmLo,d.cmHi));
      hd->Draw("E");
      struct K{int w;int col;int sty;};
      // 0 total 1 gs 2 5/2+ 3 bkg 4 P1 5 P2 6 singles
      for(int kk=0;kk<7;kk++){
        int cols[7]={kAzure+2,kViolet+1,kCyan+2,kGray+2,kGreen+3,kOrange+7,kMagenta+2};
        int stys[7]={1,1,1,7,1,1,1};
        auto* g=new TGraph();
        for(double x=B1::exLo;x<=B1::exHi;x+=0.02){
          double Xt=x-d.gsMu, v=0;
          const double* a=bp.data()+3+f.nsh;
          if(kk==0) v=B1::model(x,bp.data());
          else if(kk==1) v=std::max(0.0,bp[0])*shapeGT(x,d.gsMu,d.sgL,d.sgR,d.alL,d.alR);
          else if(kk==2) v=std::max(0.0,bp[1])*shapeGT(x,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR);
          else if(kk==3) v=std::max(0.0,bp[2])*B1::shp(Xt,bp.data()+3)*(B1::ACC?std::max(0.0,B1::ACC->Eval(Xt)):1.0);
          else { int grp=kk-4; for(int r=0;r<RK::NR;r++) if(RK::GRP[r]==grp&&d.in[r])
                   v+=std::max(0.0,a[grp])*((RK::RIDX[r]>=0)?RK::FIXRAT[RK::RIDX[r]]:1.0)*d.T[r].Eval(x);
                 if(grp==2){ for(int r=0;r<RK::NR;r++) if(RK::GRP[r]==3&&d.in[r])
                   v+=std::max(0.0,a[3])*d.T[r].Eval(x); } }
          if(v>0) g->SetPoint(g->GetN(),x,v); }
        if(g->GetN()<2) continue;
        g->SetLineColor(cols[kk]); g->SetLineStyle(stys[kk]); g->SetLineWidth(kk==0?3:2); g->Draw("L"); }
      hd->Draw("E SAME");
      auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.030);
      tx->DrawLatex(0.45,0.86,Form("bkg = p_{3}*q #upoint (E_{x}-S_{n})^{a} e^{-k(E_{x}-S_{n})}"));
      tx->DrawLatex(0.45,0.81,Form("a = %.2f    k = %.2f",bp[3],bp[4]));
      tx->DrawLatex(0.45,0.76,Form("#chi^{2}/ndf = %.2f   over %.1f to %.1f",
                                   s2/std::max(1.0,nd-(double)npar),B1::exLo,B1::exHi));
      { auto* l=new TLine(PS::SN+d.gsMu,0.4,PS::SN+d.gsMu,hd->GetMaximum()*0.3);
        l->SetLineStyle(3); l->SetLineColor(kGray+2); l->Draw();
        auto* t2=new TLatex(PS::SN+d.gsMu+0.1,hd->GetMaximum()*0.33,"S_{n}");
        t2->SetTextSize(0.030); t2->SetTextColor(kGray+2); t2->Draw(); }
      c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/bin1015_fullrange.png");
      printf("   wrote bin1015_fullrange.png\n"); }
    delete m; }
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RK::PAIRS=false; BB_QUIET=true; binbuild();
}
