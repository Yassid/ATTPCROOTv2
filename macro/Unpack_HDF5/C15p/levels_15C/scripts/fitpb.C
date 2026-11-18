#include "resfitD.C"
// THE FITS, with a background shape free PER ANGULAR BIN (states fixed).
// Each angle is an independent 6-parameter fit once the two pair ratios are fixed from the global
// solution: continuum normalisation + 4 group amplitudes + the continuum shape.  n >= 0, since a
// form factor rising with Ex is not physical.
namespace PB { int BIN=0; }
static double QS[16][6];                       // the converged parameters, per bin
static double YLD[2][2][16][2];                // [pass][code][bin][P1,P2]
static int    GPASS=0;
double pbF(const double* q){
  RK::Bin& d=RK::B[PB::BIN]; double s=0;
  for(int i=1;i<=d.h->GetNbinsX();i++){
    double x=d.h->GetBinCenter(i);
    if(x<RK::exLo||x>RK::exHi) continue;
    double n=d.h->GetBinContent(i), Xt=x-d.gsMu;
    double m=d.gsAmp*shapeGT(x,d.gsMu,d.sgL,d.sgR,d.alL,d.alR)
            +d.e5Amp*shapeGT(x,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR)
            +std::max(0.0,q[0])*RK::bshape(Xt,std::max(0.0,q[5]))
                               *(d.ga?std::max(0.0,d.ga->Eval(Xt)):1.0);
    for(int r=0;r<RK::NR;r++) if(RK::GRP[r]>=0&&d.in[r])
      m += std::max(0.0,q[1+RK::GRP[r]])*((RK::RIDX[r]>=0)?RK::FIXRAT[RK::RIDX[r]]:1.0)*d.T[r].Eval(x);
    if(m<1e-9) m=1e-9;
    s += 2.0*(m-n+(n>0? n*std::log(n/m):0.0)); }
  return s;
}
double compPB(int b,double X,int which){       // 0 tot 1 cont 2 P1 3 P2 4 gs 5 5/2+ 6 singles
  RK::Bin& d=RK::B[b]; const double* q=QS[b];
  double Xt=X-d.gsMu;
  double gs=d.gsAmp*shapeGT(X,d.gsMu,d.sgL,d.sgR,d.alL,d.alR);
  double e5=d.e5Amp*shapeGT(X,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR);
  double bg=std::max(0.0,q[0])*RK::bshape(Xt,std::max(0.0,q[5]))*(d.ga?std::max(0.0,d.ga->Eval(Xt)):1.0);
  double a[4]={0,0,0,0};
  for(int r=0;r<RK::NR;r++) if(RK::GRP[r]>=0&&d.in[r])
    a[RK::GRP[r]] += std::max(0.0,q[1+RK::GRP[r]])*((RK::RIDX[r]>=0)?RK::FIXRAT[RK::RIDX[r]]:1.0)*d.T[r].Eval(X);
  switch(which){ case 1:return bg; case 2:return a[0]; case 3:return a[1];
                 case 4:return gs; case 5:return e5; case 6:return a[2]+a[3]; }
  return gs+e5+bg+a[0]+a[1]+a[2]+a[3];
}
void fitpb(int pass=0);
void angpb();
void fitpb(int pass){
  setCM(10,5,60);
  RK::exLo=3.5; RK::exHi=9.0; RK::ONETERM=true; RK::BKG=RK::HALON;
  RK::TIE12=false; RK::PAIRS=true; RK::NPERBIN=false; RK::DELTA=0.0;
  for(int i=0;i<RK::NR;i++) RK::USE[i]=(i<6);
  Cfg b0=CFG[0], b1=CFG[1];
  GPASS=pass;
  gStyle->SetOptStat(0);
  TCanvas* c=new TCanvas(Form("c%d",pass),"c",2000,1500); c->Divide(5,4,0.001,0.001);
  int col[2]={kAzure+2,kRed+1};
  for(int ic=0;ic<2;ic++){
    CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; if(ic==0) CFG[ic].chi2=5.0;
    VZLO = pass? 350 : 300; BB_QUIET=true; binbuild();
    RK::SEED.clear();
    rkRun(CFG[ic],NCM,false);
    double nglob=RK::LASTP0;
    for(int r=0;r<RK::NRAT;r++) RK::FIXRAT[r]=RK::LASTP[1+NCM*RK::nper()+r];
    printf("\n================ %s : per-angle background, states fixed ================\n",CFG[ic].nm);
    printf("   pair ratios %.2f / %.2f (global)\n",RK::FIXRAT[0],RK::FIXRAT[1]);
    printf("   %-10s %8s %8s %10s %10s %10s\n","theta_CM","chi2/ndf","n","P1 yield","P2 yield","cont@6");
    double mc=0; int nb=0;
    for(int b=0;b<NCM;b++){
      RK::Bin& d=RK::B[b]; if(!d.h) continue;
      PB::BIN=b; RK::ONLYBIN=b;
      int o=1+b*RK::nper();
      double q0[6]={RK::LASTP[o+1],0,0,0,0,nglob};
      for(int g=0;g<RK::NGRP;g++) q0[1+g]=RK::LASTP[o+2+g];
      auto* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
      m->SetMaxFunctionCalls(200000); m->SetTolerance(0.01); m->SetPrintLevel(0); m->SetStrategy(2);
      ROOT::Math::Functor fn(&pbF,6); m->SetFunction(fn);
      m->SetLimitedVariable(0,"Bg",q0[0],std::max(1e-6,q0[0]*0.01),0,std::max(1.0,q0[0]*1e3));
      for(int g=0;g<4;g++)
        m->SetLimitedVariable(1+g,Form("A%d",g),q0[1+g],std::max(1e-6,q0[1+g]*0.01+1e-3),
                              0,std::max(1.0,q0[1+g]*1e3+1e3));
      m->SetLimitedVariable(5,"n",nglob,0.02,0.0,9.0);
      m->Minimize();
      for(int j=0;j<6;j++) QS[b][j]=m->X()[j];
      delete m; RK::ONLYBIN=-1;
      // draw
      c->cd(1+ic*10+b); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.14); gPad->SetLogy();
      auto* hd=(TH1D*)d.h->Clone(Form("p%d%d",ic,b)); hd->SetDirectory(nullptr);
      hd->GetXaxis()->SetRangeUser(-2.0,11.0);
      double mx=0; for(int i=1;i<=hd->GetNbinsX();i++){ double x=hd->GetBinCenter(i);
        if(x>=-2&&x<=11) mx=std::max(mx,hd->GetBinContent(i)); }
      hd->SetMaximum(6.0*std::max(1.0,mx)); hd->SetMinimum(0.35);
      hd->SetMarkerStyle(20); hd->SetMarkerSize(0.45); hd->SetLineColor(kBlack);
      hd->SetTitle(Form("%s  #theta_{CM} %.0f-%.0f#circ;E_{x} [MeV];counts / bin",
                        CFG[ic].nm,d.cmLo,d.cmHi));
      hd->GetXaxis()->SetTitleSize(0.055); hd->GetYaxis()->SetTitleSize(0.055);
      hd->GetXaxis()->SetLabelSize(0.05);  hd->GetYaxis()->SetLabelSize(0.05);
      hd->Draw("E");
      { auto* bx=new TBox(-2.0,0.35,RK::exLo,hd->GetMaximum()); bx->SetFillColorAlpha(kGray+1,0.22);
        bx->SetLineColor(0); bx->Draw(); hd->Draw("E SAME"); }
      struct K{int w;int col;int sty;int wid;bool full;};
      K ks[6]={{4,kViolet+1,1,2,true},{5,kCyan+2,1,2,true},{1,kGray+2,7,2,false},
               {2,kGreen+3,1,2,false},{3,kOrange+7,1,2,false},{0,col[ic],1,3,false}};
      for(auto&kk:ks){ auto* g=new TGraph();
        double x0=kk.full?-2.0:RK::exLo, x1=kk.full?11.0:RK::exHi;
        for(double x=x0;x<=x1+1e-9;x+=0.02){ double v=compPB(b,x,kk.w); if(v>0) g->SetPoint(g->GetN(),x,v); }
        if(g->GetN()<2) continue;
        g->SetLineColor(kk.col); g->SetLineStyle(kk.sty); g->SetLineWidth(kk.wid); g->Draw("L"); }
      hd->Draw("E SAME");
      double s=0; int nd=0;
      for(int i=1;i<=hd->GetNbinsX();i++){ double x=hd->GetBinCenter(i),n=hd->GetBinContent(i);
        if(x<RK::exLo||x>RK::exHi) continue;
        double mm=compPB(b,x,0); if(mm<1e-9) mm=1e-9; s+=(n-mm)*(n-mm)/std::max(mm,1.0); nd++; }
      double ndf=std::max(1.0,nd-6.0);
      auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.058);
      tx->DrawLatex(0.52,0.87,Form("#chi^{2}/ndf %.2f",s/ndf));
      tx->SetTextSize(0.044); tx->SetTextColor(kGray+2);
      tx->DrawLatex(0.52,0.805,Form("n = %.2f",QS[b][5]));
      tx->DrawLatex(0.52,0.745,Form("N %.0f",d.N));
      double p1=compPB(b,(RK::EPP6[0]-RK::PPGS6)+d.gsMu,2);
      double p2=compPB(b,(RK::EPP6[3]-RK::PPGS6)+d.gsMu,3);
      YLD[GPASS][ic][b][0]=p1; YLD[GPASS][ic][b][1]=p2;
      printf("   %3.0f-%-6.0f %8.2f %8.2f %10.1f %10.1f %10.1f\n",
             d.cmLo,d.cmHi,s/ndf,QS[b][5],p1,p2,compPB(b,6.0,1));
      mc+=s/ndf; nb++; }
    printf("   mean chi2/ndf %.2f\n",nb?mc/nb:0);
  }
  if(!pass){ c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/fit_perangle_bkg.png");
    printf("\n   wrote fit_perangle_bkg.png\n"); }
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RK::PAIRS=false; BB_QUIET=true; binbuild();
}

// The angular distributions that come out of the per-angle-background fit, with the stability test
// that has decided every earlier question: refit under a deliberately small change (VZLO 300->350,
// about 5 % of the sample) and show the spread as part of the error.
void angpb(){
  fitpb(0); fitpb(1);
  gStyle->SetOptStat(0);
  TCanvas* c=new TCanvas("ca","ca",1150,540); c->Divide(2,1);
  const char* pn[2]={"P1  (4.43 + 4.93)","P2  (5.98 + 6.68)"};
  int col[2]={kAzure+2,kRed+1}, mk[2]={20,21};
  for(int ip=0;ip<2;ip++){
    c->cd(ip+1); gPad->SetLeftMargin(0.15); gPad->SetBottomMargin(0.14);
    gPad->SetGridy(); gPad->SetLogy();
    double ymax=0,ymin=1e30;
    TGraphErrors* g[2];
    for(int ic=0;ic<2;ic++){
      g[ic]=new TGraphErrors();
      for(int b=0;b<NCM;b++){
        double y0=YLD[0][ic][b][ip], y1=YLD[1][ic][b][ip];
        if(y0<=0.3) continue;                       // below one count: not a measurement
        double sr=std::fabs(y1/0.94-y0);
        double e=std::max(std::sqrt(std::max(1.0,y0)),sr);
        int p=g[ic]->GetN();
        g[ic]->SetPoint(p,TCM0+(b+0.5)*TCMSTEP,y0); g[ic]->SetPointError(p,TCMSTEP/2,e);
        ymax=std::max(ymax,y0+e); ymin=std::min(ymin,std::max(0.3,y0-e)); }
      g[ic]->SetMarkerStyle(mk[ic]); g[ic]->SetMarkerSize(1.3);
      g[ic]->SetMarkerColor(col[ic]); g[ic]->SetLineColor(col[ic]); g[ic]->SetLineWidth(2); }
    auto* fr=gPad->DrawFrame(TCM0-2,0.4*ymin,TCM0+NCM*TCMSTEP+2,4.0*ymax);
    fr->SetTitle(Form("%s;#theta_{CM} [deg]  (projectile convention);yield [counts / bin]",pn[ip]));
    for(int ic=0;ic<2;ic++) g[ic]->Draw("P");
    auto* lg=new TLegend(0.60,0.74,0.92,0.88); lg->SetBorderSize(0); lg->SetFillStyle(0);
    lg->AddEntry(g[0],"ATTPCROOT","pl"); lg->AddEntry(g[1],"Spyral","pl"); lg->Draw();
    auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.028); tx->SetTextColor(kGray+3);
    tx->DrawLatex(0.18,0.21,"#splitline{v_{z}#in[300,900], #chi^{2}<5, fit 3.5-9,}{background shape free PER ANGLE}");
    tx->DrawLatex(0.18,0.145,"errors include the refit spread");
  }
  c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/angdist_perangle_bkg.png");
  printf("\n   wrote angdist_perangle_bkg.png\n");
}
