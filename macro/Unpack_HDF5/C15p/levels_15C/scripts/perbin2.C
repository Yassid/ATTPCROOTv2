#include "resfitD.C"
// ONE ANGULAR BIN AT A TIME.  States fixed; the continuum shape free per bin.
// With the two pair ratios held at their global-fit values each angle is an independent
// 6-parameter problem (continuum normalisation + 4 group amplitudes + the shape), which converges
// where the 73-parameter simultaneous fit did not.
namespace PB { int BIN=0; TH1D* H=nullptr; }
double pbFcn(const double* q){
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
      m += std::max(0.0,q[1+RK::GRP[r]]) * ((RK::RIDX[r]>=0)?RK::FIXRAT[RK::RIDX[r]]:1.0)
           * d.T[r].Eval(x);
    if(m<1e-9) m=1e-9;
    s += 2.0*(m-n+(n>0? n*std::log(n/m):0.0)); }
  return s;
}
void perbin2(){
  setCM(10,5,60);
  RK::exLo=3.5; RK::exHi=9.0; RK::ONETERM=true; RK::BKG=RK::HALON;
  RK::TIE12=false; RK::PAIRS=true; RK::NPERBIN=false; RK::DELTA=0.0;
  for(int i=0;i<RK::NR;i++) RK::USE[i]=(i<6);
  Cfg b0=CFG[0], b1=CFG[1];
  for(int ic=0;ic<2;ic++){
    CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; if(ic==0) CFG[ic].chi2=5.0;
    VZLO=300; BB_QUIET=true; binbuild();
    RK::SEED.clear();
    double Lg=rkRun(CFG[ic],NCM,false);          // global shape, to set the ratios
    double nglob=RK::LASTP0;
    for(int r=0;r<RK::NRAT;r++) RK::FIXRAT[r]=RK::LASTP[1+NCM*RK::nper()+r];
    printf("\n================ %s : background shape per ANGLE ================\n",CFG[ic].nm);
    printf("   global shape fit: -2lnL %.1f   n = %.2f   ratios %.2f / %.2f\n",
           Lg,nglob,RK::FIXRAT[0],RK::FIXRAT[1]);
    printf("   %-10s %9s %9s %9s %9s %9s\n",
           "theta_CM","chi2 glob","chi2 /bin","n fitted","-2lnL g","-2lnL b");
    double mg=0,mb=0; int nb=0; double sLg=0,sLb=0;
    for(int b=0;b<NCM;b++){
      RK::Bin& d=RK::B[b]; if(!d.h) continue;
      PB::BIN=b; RK::ONLYBIN=b;
      auto chib=[&](const double* q,double shp){
        double s=0; int nd=0;
        for(int i=1;i<=d.h->GetNbinsX();i++){
          double x=d.h->GetBinCenter(i),n=d.h->GetBinContent(i);
          if(x<RK::exLo||x>RK::exHi) continue;
          double Xt=x-d.gsMu;
          double m=d.gsAmp*shapeGT(x,d.gsMu,d.sgL,d.sgR,d.alL,d.alR)
                  +d.e5Amp*shapeGT(x,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR)
                  +std::max(0.0,q[0])*RK::bshape(Xt,shp)*(d.ga?std::max(0.0,d.ga->Eval(Xt)):1.0);
          for(int r=0;r<RK::NR;r++) if(RK::GRP[r]>=0&&d.in[r])
            m += std::max(0.0,q[1+RK::GRP[r]])*((RK::RIDX[r]>=0)?RK::FIXRAT[RK::RIDX[r]]:1.0)
                 *d.T[r].Eval(x);
          if(m<1e-9) m=1e-9; s+=(n-m)*(n-m)/std::max(m,1.0); nd++; }
        return s/std::max(1.0,nd-6.0); };
      // seed from the global solution for this bin
      int o=1+b*RK::nper();
      double q0[6]={RK::LASTP[o+1],0,0,0,0,nglob};
      for(int g=0;g<RK::NGRP;g++) q0[1+g]=RK::LASTP[o+2+g];
      double cg=chib(q0,nglob);
      double Lgb=pbFcn(q0);
      auto* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
      m->SetMaxFunctionCalls(200000); m->SetTolerance(0.01); m->SetPrintLevel(0); m->SetStrategy(2);
      ROOT::Math::Functor fn(&pbFcn,6); m->SetFunction(fn);
      m->SetLimitedVariable(0,"Bg",q0[0],std::max(1e-6,q0[0]*0.01),0,std::max(1.0,q0[0]*1e3));
      for(int g=0;g<4;g++)
        m->SetLimitedVariable(1+g,Form("A%d",g),q0[1+g],std::max(1e-6,q0[1+g]*0.01+1e-3),
                              0,std::max(1.0,q0[1+g]*1e3+1e3));
      m->SetLimitedVariable(5,"n",nglob,0.02,-4.0,9.0);   // allow NEGATIVE: 6-7 bins were pinned at the 0 bound
      m->Minimize();
      const double* q=m->X();
      double cb=chib(q,q[5]); double Lb=m->MinValue();
      printf("   %3.0f-%-6.0f %9.2f %9.2f %9.2f %9.1f %9.1f%s\n",
             d.cmLo,d.cmHi,cg,cb,q[5],Lgb,Lb,(Lb>Lgb+0.01)?"  <- WORSE, not converged":"");
      mg+=cg; mb+=cb; nb++; sLg+=Lgb; sLb+=Lb;
      delete m; RK::ONLYBIN=-1; }
    printf("   mean chi2  %.2f (one shape)  ->  %.2f (per angle)   total -2lnL %.1f -> %.1f  gain %+.1f for %d par\n",
           nb?mg/nb:0, nb?mb/nb:0, sLg, sLb, sLg-sLb, NCM);
  }
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RK::PAIRS=false; BB_QUIET=true; binbuild();
}
