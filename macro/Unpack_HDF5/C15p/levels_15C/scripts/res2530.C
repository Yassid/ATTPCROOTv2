#include "fitpb.C"
// 25-30 is the worst bin in BOTH codes and it is where the continuum exponent transitions from
// ~1 to 0.  Is the power law too rigid there, or is something else wrong?  Print the pulls.
void res2530(){
  fitpb(0);
  for(int ic=0;ic<2;ic++){
    // refit that one bin so QS/RK::B belong to this code
    Cfg b0=CFG[0], b1=CFG[1];
    CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; if(ic==0) CFG[ic].chi2=5.0;
    VZLO=300; BB_QUIET=true; binbuild();
    RK::SEED.clear(); rkRun(CFG[ic],NCM,false);
    double nglob=RK::LASTP0;
    for(int r=0;r<RK::NRAT;r++) RK::FIXRAT[r]=RK::LASTP[1+NCM*RK::nper()+r];
    int b=3; PB::BIN=b; RK::ONLYBIN=b;
    RK::Bin& d=RK::B[b];
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
    printf("\n======== %s  theta_CM 25-30  (n = %.2f) ========\n",CFG[ic].nm,QS[b][5]);
    printf("   %6s %8s %9s %8s   %8s %8s %8s  %s\n",
           "Ex","data","model","pull","cont","P1","P2","");
    for(int i=1;i<=d.h->GetNbinsX();i++){
      double x=d.h->GetBinCenter(i);
      if(x<RK::exLo||x>RK::exHi) continue;
      double n=d.h->GetBinContent(i), mm=compPB(b,x,0);
      if(mm<1e-9) mm=1e-9;
      double pull=(n-mm)/std::sqrt(std::max(mm,1.0));
      char bar[40]; int k=0; int np=(int)std::min(16.0,std::fabs(pull)*3);
      for(int j=0;j<np;j++) bar[k++]= pull>0?'+':'-';
      bar[k]=0;
      printf("   %6.2f %8.0f %9.1f %8.2f   %8.1f %8.1f %8.1f  %s\n",
             x,n,mm,pull,compPB(b,x,1),compPB(b,x,2),compPB(b,x,3),bar); }
    CFG[0]=b0; CFG[1]=b1;
  }
  VZLO=0; RK::PAIRS=false; BB_QUIET=true; binbuild();
}
