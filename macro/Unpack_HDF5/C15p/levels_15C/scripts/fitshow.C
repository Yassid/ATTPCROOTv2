#include "angfinal.C"
// The fits as they now stand: peak EXCLUDED from the fit (so no g.s./5/2+ parameters -- 7 instead
// of 11, which is what finally made them converge), background free per angle and defined from
// Sn, states fixed.  Data drawn over the full range; the model only where it is fitted.
static double QQ[2][16][12]; static int QST[2][16]; static double QE[2][16], QX[2][16];
void fitshow(){
  setCM(10,5,60);
  RK::exLo=3.5; RK::exHi=9.0; RK::ONETERM=true; RK::BKG=RK::HALON;
  RK::TIE12=false; RK::PAIRS=true; RK::DELTA=0.0; B1::exLo=-1.5;
  for(int i=0;i<RK::NR;i++) RK::USE[i]=(i<6);
  Cfg b0=CFG[0], b1=CFG[1];
  gStyle->SetOptStat(0);
  TCanvas* c=new TCanvas("c","c",2000,1500); c->Divide(5,4,0.001,0.001);
  int col[2]={kAzure+2,kRed+1};
  for(int ic=0;ic<2;ic++){
    CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; if(ic==0) CFG[ic].chi2=5.0;
    VZLO=300; BB_QUIET=true; binbuild();
    RK::SEED.clear(); rkRun(CFG[ic],NCM,false);
    for(int r=0;r<RK::NRAT;r++) RK::FIXRAT[r]=RK::LASTP[1+NCM*RK::nper()+r];
    B1::EB=CFG[ic].eb;
    for(int b=0;b<NCM;b++){
      RK::Bin& d=RK::B[b]; if(!d.h) continue;
      B1::BIN=b; RK::ONLYBIN=b; B1::PEAKFREE=true;
      Pars A=parsOf(CFG[ic]); A.vzLo=300; Win w{14,86,CFG[ic].keLo,CFG[ic].keHi};
      B1::ACC=new TGraph();
      for(int i=0;i<=460;i++){ double X=i*0.05;
        B1::ACC->SetPoint(i,X,accBin(A,CFG[ic].eb,X,d.cmLo,d.cmHi,w)); }
      const int npar=11; B1::NPAR=npar; B1::NSH=2; B1::FORM=3;
      std::vector<double> best(npar,0); double bestL=1e18; int st=-1; double edm=-1;
      const double AG[9]={0.0,0.25,0.5,0.75,1.0,1.25,1.5,2.0,3.0};
      for(double av:AG){
        auto* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
        m->SetMaxFunctionCalls(2000000); m->SetTolerance(0.001);
        m->SetPrintLevel(0); m->SetStrategy(2);
        ROOT::Math::Functor fn(&B1::fcn,npar); m->SetFunction(fn);
        double pk=d.h->GetBinContent(d.h->GetMaximumBin());
        m->SetLimitedVariable(0,"Ags",pk,pk*0.01,0,pk*10);
        m->SetLimitedVariable(1,"A52",0.1*pk,pk*0.001,0,pk*5);
        m->SetLimitedVariable(2,"B",50,1,0,1e7);
        m->SetFixedVariable(3,"a",av);
        m->SetLimitedVariable(4,"k",0.5,0.02,0.0,4.0);
        for(int g=0;g<4;g++) m->SetLimitedVariable(5+g,Form("s%d",g),4.0,0.05,-80,80);
        m->SetLimitedVariable(npar-2,"mu",d.gsMu,0.005,d.gsMu-0.6,d.gsMu+0.6);
        m->SetLimitedVariable(npar-1,"sig",0.4,0.005,0.08,1.5);
        m->Minimize();
        for(int t=0;t<3&&m->Status()!=0;t++) m->Minimize();
        m->Hesse();
        if(m->MinValue()<bestL){ bestL=m->MinValue(); best.assign(m->X(),m->X()+npar);
          st=m->Status(); edm=m->Edm(); }
        delete m; }
      { auto* mf=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
        mf->SetMaxFunctionCalls(2000000); mf->SetTolerance(0.001);
        mf->SetPrintLevel(0); mf->SetStrategy(2);
        ROOT::Math::Functor fnf(&B1::fcn,npar); mf->SetFunction(fnf);
        double pk=d.h->GetBinContent(d.h->GetMaximumBin());
        mf->SetLimitedVariable(0,"Ags",best[0],std::max(1e-3,best[0]*0.01),0,pk*10);
        mf->SetLimitedVariable(1,"A52",best[1],std::max(1e-3,best[1]*0.01+1e-3),0,pk*5);
        mf->SetLimitedVariable(2,"B",best[2],std::max(1e-3,best[2]*0.01),0,1e7);
        mf->SetLimitedVariable(3,"a",std::max(0.02,best[3]),0.02,0.0,6.0);
        mf->SetLimitedVariable(4,"k",best[4],0.02,0.0,4.0);
        for(int g=0;g<4;g++) mf->SetLimitedVariable(5+g,Form("s%d",g),best[5+g],
                                 std::max(1e-3,std::fabs(best[5+g])*0.01+1e-3),-80,80);
        mf->SetLimitedVariable(npar-2,"mu",best[npar-2],0.005,d.gsMu-0.6,d.gsMu+0.6);
        mf->SetLimitedVariable(npar-1,"sig",best[npar-1],0.005,0.08,1.5);
        mf->Minimize();
        for(int t=0;t<3&&mf->Status()!=0;t++) mf->Minimize();
        mf->Hesse();
        if(mf->MinValue()<bestL){ bestL=mf->MinValue(); best.assign(mf->X(),mf->X()+npar);
          st=mf->Status(); edm=mf->Edm(); }
        delete mf; }
      for(int j=0;j<npar;j++) QQ[ic][b][j]=best[j];
      QST[ic][b]=st; QE[ic][b]=edm;
      const double* p=best.data();
      double s2=0; int nd=0;
      for(int i=1;i<=d.h->GetNbinsX();i++){ double x=d.h->GetBinCenter(i),n=d.h->GetBinContent(i);
        if(x<B1::exLo||x>B1::exHi) continue;
        double mm=B1::model(x,p); if(mm<1e-9) mm=1e-9; s2+=(n-mm)*(n-mm)/std::max(mm,1.0); nd++; }
      QX[ic][b]=s2/std::max(1.0,nd-7.0);
      // ---- draw
      c->cd(1+ic*10+b); gPad->SetLeftMargin(0.13); gPad->SetBottomMargin(0.14); gPad->SetLogy();
      auto* hd=(TH1D*)d.h->Clone(Form("f%d%d",ic,b)); hd->SetDirectory(nullptr);
      hd->GetXaxis()->SetRangeUser(-2.0,11.0);
      double mx=0; for(int i=1;i<=hd->GetNbinsX();i++){ double x=hd->GetBinCenter(i);
        if(x>=-2&&x<=11) mx=std::max(mx,hd->GetBinContent(i)); }
      hd->SetMaximum(6*std::max(1.0,mx)); hd->SetMinimum(0.35);
      hd->SetMarkerStyle(20); hd->SetMarkerSize(0.45); hd->SetLineColor(kBlack);
      hd->SetTitle(Form("%s  #theta_{CM} %.0f-%.0f#circ;E_{x} [MeV];counts / 100 keV",
                        CFG[ic].nm,d.cmLo,d.cmHi));
      hd->GetXaxis()->SetTitleSize(0.055); hd->GetYaxis()->SetTitleSize(0.055);
      hd->GetXaxis()->SetLabelSize(0.05);  hd->GetYaxis()->SetLabelSize(0.05);
      hd->Draw("E");

      // components: background, P1, P2, singles, total
      for(int kk=0;kk<7;kk++){
        int cols[7]={kGray+2,kGreen+3,kOrange+7,kMagenta+2,col[ic],kViolet+1,kCyan+2};
        int stys[7]={7,1,1,1,1,1,1};
        auto* g=new TGraph();
        for(double x=B1::exLo;x<=B1::exHi+1e-9;x+=0.02){
          double mu=p[npar-2], sg=std::max(0.05,p[npar-1]);
          double Xt=x-mu, v=0; const double* a2=p+5;
          if(kk==4) v=B1::model(x,p);
          else if(kk==5){ double t=(x-mu)/sg; v=std::max(0.0,p[0])*std::exp(-0.5*t*t); }
          else if(kk==6){ double t=(x-mu-0.740)/sg; v=std::max(0.0,p[1])*std::exp(-0.5*t*t); }
          else if(kk==0) v=std::max(0.0,p[2])*B1::shp(Xt,p+3)*(B1::ACC?std::max(0.0,B1::ACC->Eval(Xt)):1.0);
          else { int grp=kk-1;
            for(int r=0;r<RK::NR;r++) if(RK::GRP[r]==grp&&d.in[r])
              v+=a2[grp]*a2[grp]*((RK::RIDX[r]>=0)?RK::FIXRAT[RK::RIDX[r]]:1.0)*d.T[r].Eval(x);
            if(grp==2) for(int r=0;r<RK::NR;r++) if(RK::GRP[r]==3&&d.in[r])
              v+=a2[3]*a2[3]*d.T[r].Eval(x); }
          if(v>0) g->SetPoint(g->GetN(),x,v); }
        if(g->GetN()<2) continue;
        g->SetLineColor(cols[kk]); g->SetLineStyle(stys[kk]); g->SetLineWidth(kk==4?3:2); g->Draw("L"); }
      hd->Draw("E SAME");
      auto* tx=new TLatex(); tx->SetNDC(); tx->SetTextSize(0.057);
      tx->DrawLatex(0.50,0.88,Form("#chi^{2}/ndf %.2f",QX[ic][b]));
      tx->SetTextSize(0.044);
      tx->SetTextColor(st==0?kGreen+3:kRed+1);
      tx->DrawLatex(0.50,0.815,Form("status %d",st));
      tx->SetTextColor(kGray+2);
      tx->DrawLatex(0.50,0.755,Form("EDM %.0e",edm));
      tx->DrawLatex(0.50,0.695,Form("a %.2f  k %.2f",p[3],p[4]));
      RK::ONLYBIN=-1; }
  }
  c->SaveAs("/home/yassid/fair_install/ATTPCROOTv2-OpenKF/macro/Unpack_HDF5/C15p/levels_15C/figures/fits_fullrange.png");
  printf("   wrote fits_fullrange.png\n");
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RK::PAIRS=false; BB_QUIET=true; binbuild();
}
