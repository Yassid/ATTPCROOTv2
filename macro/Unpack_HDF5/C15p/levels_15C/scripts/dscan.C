#include "resfitD.C"
// Is our Ex scale simply OFFSET?  Shift every state by one common DELTA and refit everything.
// A calibration offset moves all states together; per-state tuning cannot repair it.
void dscan(){
  setCM(10,5,60);
  RK::exLo=3.5; RK::exHi=9.0; RK::ONETERM=true; RK::BKG=RK::HALON;
  RK::TIE12=false; RK::PAIRS=true; RK::NPERBIN=false;
  for(int i=0;i<RK::NR;i++) RK::USE[i]=(i<6);
  Cfg b0=CFG[0], b1=CFG[1];
  for(int ic=0;ic<2;ic++){
    CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; if(ic==0) CFG[ic].chi2=5.0;
    VZLO=300; BB_QUIET=true; binbuild();
    printf("\n================ %s : one global energy shift ================\n",CFG[ic].nm);
    printf("   %8s %11s %10s %10s\n","DELTA","-2lnL","gain","mean chi2");
    double L0=0,best=1e18,bD=0;
    for(double D=-0.40;D<=0.401;D+=0.05){
      RK::DELTA=D;
      double L=rkRun(CFG[ic],NCM,false);
      double mc=0; int nb=0;
      for(int b=0;b<NCM;b++){ RK::Bin& d=RK::B[b]; if(!d.h) continue;
        double s=0; int nd=0;
        for(int i=1;i<=d.h->GetNbinsX();i++){
          double x=d.h->GetBinCenter(i),n=d.h->GetBinContent(i);
          if(x<RK::exLo||x>RK::exHi) continue;
          double Xt=x-d.gsMu;
          double m=d.gsAmp*shapeGT(x,d.gsMu,d.sgL,d.sgR,d.alL,d.alR)
                  +d.e5Amp*shapeGT(x,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR)
                  +RK::LASTBG[b]*RK::bshape(Xt,RK::LASTP0)*(d.ga?std::max(0.0,d.ga->Eval(Xt)):1.0);
          for(int r=0;r<RK::NR;r++) if(RK::GRP[r]>=0&&d.in[r]) m+=RK::LASTA[b][r]*d.T[r].Eval(x);
          if(m<1e-9) m=1e-9; s+=(n-m)*(n-m)/std::max(m,1.0); nd++; }
        mc+=s/std::max(1.0,nd-5.3); nb++; }
      if(std::fabs(D)<1e-9) L0=L;
      if(L<best){best=L;bD=D;}
      printf("   %+8.2f %11.1f %10s %10.2f%s\n",D,L,
             (L0>0?Form("%+.1f",L0-L):"-"),nb?mc/nb:0,(std::fabs(D)<1e-9)?"   <- current":""); }
    printf("   best shift %+.2f MeV, gain %.1f over no shift\n",bD,L0-best);
  }
  RK::DELTA=0; CFG[0]=b0; CFG[1]=b1; VZLO=0; RK::PAIRS=false; BB_QUIET=true; binbuild();
}
