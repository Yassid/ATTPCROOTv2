#include "resfitD.C"
// Bin 15-20 is the richest and the worst fitted in BOTH codes (chi2/ndf 5.5 and 6.7), under every
// resonance scheme tried.  WHERE in Ex does the model fail?  Print the pull per 100 keV bin.
void resid1520(){
  setCM(10,5,60);
  RK::exLo=3.0; RK::exHi=9.0; RK::ONETERM=true; RK::BKG=RK::HALON;
  RK::TIE12=false; RK::PAIRS=true;
  bool keep[6]={true,true,true,true,false,false};
  for(int i=0;i<RK::NR;i++) RK::USE[i]=keep[i];
  Cfg b0=CFG[0], b1=CFG[1];
  for(int ic=0;ic<2;ic++){
    CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; VZLO=300; BB_QUIET=true; binbuild();
    rkRun(CFG[ic],NCM,false);
    int b=1;  // theta_CM 15-20
    RK::Bin& d=RK::B[b];
    printf("\n======== %s  theta_CM 15-20 : residuals ========\n",CFG[ic].nm);
    printf("   %6s %8s %9s %8s   %s\n","Ex","data","model","pull","");
    double worst=0; double wx=0;
    for(int i=1;i<=d.h->GetNbinsX();i++){
      double x=d.h->GetBinCenter(i), n=d.h->GetBinContent(i);
      if(x<RK::exLo||x>RK::exHi) continue;
      double Xt=x-d.gsMu;
      double m=d.gsAmp*shapeGT(x,d.gsMu,d.sgL,d.sgR,d.alL,d.alR)
              +d.e5Amp*shapeGT(x,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR)
              +RK::LASTBG[b]*RK::bshape(Xt,RK::LASTP0)*(d.ga?std::max(0.0,d.ga->Eval(Xt)):1.0);
      for(int r=0;r<RK::NR;r++) if(RK::GRP[r]>=0&&d.in[r]) m+=RK::LASTA[b][r]*d.T[r].Eval(x);
      if(m<1e-9) m=1e-9;
      double pull=(n-m)/std::sqrt(std::max(1.0,m));
      if(std::fabs(pull)>std::fabs(worst)){ worst=pull; wx=x; }
      char bar[48]; int k=0; int np=(int)std::min(20.0,std::fabs(pull)*3);
      for(int j=0;j<np;j++) bar[k++]= pull>0?'+':'-';
      bar[k]=0;
      printf("   %6.2f %8.0f %9.1f %8.2f   %s\n",x,n,m,pull,bar); }
    printf("   worst pull %+.2f at Ex = %.2f\n",worst,wx);
  }
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RK::PAIRS=false;
  for(int i=0;i<RK::NR;i++) RK::USE[i]=true; BB_QUIET=true; binbuild();
}
