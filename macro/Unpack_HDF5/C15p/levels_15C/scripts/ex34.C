#include "resfitD.C"
// Scan a seventh state through Ex 3.0-4.0 and ask whether the richest bin recovers.
// Threshold for 10 free amplitudes is Delta(-2lnL) ~ 18-25.
void ex34(){
  setCM(10,5,60);
  RK::exLo=3.0; RK::exHi=9.0; RK::ONETERM=true; RK::BKG=RK::HALON;
  RK::TIE12=false; RK::PAIRS=true;
  Cfg b0=CFG[0], b1=CFG[1];
  auto chi1520=[&](int b){ RK::Bin& d=RK::B[b]; double s=0; int nd=0;
    for(int i=1;i<=d.h->GetNbinsX();i++){ double x=d.h->GetBinCenter(i),n=d.h->GetBinContent(i);
      if(x<RK::exLo||x>RK::exHi||n<5) continue;
      double Xt=x-d.gsMu;
      double m=d.gsAmp*shapeGT(x,d.gsMu,d.sgL,d.sgR,d.alL,d.alR)
              +d.e5Amp*shapeGT(x,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR)
              +RK::LASTBG[b]*RK::bshape(Xt,RK::LASTP0)*(d.ga?std::max(0.0,d.ga->Eval(Xt)):1.0);
      for(int r=0;r<RK::NR;r++) if(RK::GRP[r]>=0&&d.in[r]) m+=RK::LASTA[b][r]*d.T[r].Eval(x);
      if(m<1e-9) m=1e-9; s+=(n-m)*(n-m)/n; nd++; }
    return s/std::max(1,nd-5); };
  for(int ic=0;ic<2;ic++){
    CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; VZLO=300; BB_QUIET=true; binbuild();
    printf("\n======== %s : a seventh state below the doublet ========\n",CFG[ic].nm);
    bool off[7]={true,true,true,true,false,false,false};
    for(int i=0;i<RK::NR;i++) RK::USE[i]=off[i];
    double L0=rkRun(CFG[ic],NCM,false); double c0=chi1520(1);
    printf("   %-22s  %10s %12s %10s\n","seventh state","-2lnL","gain","chi2 15-20");
    printf("   %-22s  %10.1f %12s %10.2f\n","none",L0,"-",c0);
    double best=1e18, bE=0;
    for(double E=3.10;E<=4.01;E+=0.15){
      RK::EPP6[6]=E+RK::PPGS6;                 // EPP is on the (p,p') scale
      bool on[7]={true,true,true,true,false,false,true};
      for(int i=0;i<RK::NR;i++) RK::USE[i]=on[i];
      double L=rkRun(CFG[ic],NCM,false); double c=chi1520(1);
      printf("   Ex = %-16.2f  %10.1f %12.1f %10.2f%s\n",E,L,L0-L,c,(L0-L>25)?"  <- significant":"");
      if(L<best){best=L;bE=E;} }
    printf("   best at Ex = %.2f, gain %.1f for 10 parameters\n",bE,L0-best);
  }
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RK::PAIRS=false;
  for(int i=0;i<RK::NR;i++) RK::USE[i]=(i<6); BB_QUIET=true; binbuild();
}
