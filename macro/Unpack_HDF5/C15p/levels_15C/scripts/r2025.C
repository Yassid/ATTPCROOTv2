#include "angfinal.C"
// ATTPCROOT 20-25 is the worst surviving bin (chi2R 4.45) while Spyral's same bin is 0.86.
// It also carries the P2 peak, so the angular distribution hinges on it.  Where does it fail?
void r2025(){
  setCM(10,5,60);
  RK::exLo=3.5; RK::exHi=9.0; RK::ONETERM=true; RK::BKG=RK::HALON;
  RK::TIE12=false; RK::PAIRS=true; RK::DELTA=0.0;
  for(int i=0;i<RK::NR;i++) RK::USE[i]=(i<6);
  Cfg b0=CFG[0], b1=CFG[1];
  for(int ic=0;ic<2;ic++){
    CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; if(ic==0) CFG[ic].chi2=5.0;
    VZLO=300; BB_QUIET=true; binbuild();
    RK::SEED.clear(); rkRun(CFG[ic],NCM,false);
    for(int r=0;r<RK::NRAT;r++) RK::FIXRAT[r]=RK::LASTP[1+NCM*RK::nper()+r];
    int b=2; B1::BIN=b; B1::EB=CFG[ic].eb; RK::ONLYBIN=b;
    RK::Bin& d=RK::B[b];
    Pars A=parsOf(CFG[ic]); A.vzLo=300; Win w{14,86,CFG[ic].keLo,CFG[ic].keHi};
    B1::ACC=new TGraph();
    for(int i=0;i<=460;i++){ double X=i*0.05;
      B1::ACC->SetPoint(i,X,accBin(A,CFG[ic].eb,X,d.cmLo,d.cmHi,w)); }
    const int npar=3+2+4+2; B1::NPAR=npar; B1::NSH=2; B1::FORM=3;
    double pk=d.h->GetBinContent(d.h->GetMaximumBin());
    std::vector<double> best(npar,0); double bestL=1e18;
    const double AG[9]={0.0,0.25,0.5,0.75,1.0,1.25,1.5,2.0,3.0};
    for(double av:AG){
      auto* m=ROOT::Math::Factory::CreateMinimizer("Minuit2","Migrad");
      m->SetMaxFunctionCalls(500000); m->SetTolerance(0.01); m->SetPrintLevel(0); m->SetStrategy(2);
      ROOT::Math::Functor fn(&B1::fcn,npar); m->SetFunction(fn);
      m->SetLimitedVariable(0,"Ags",pk,pk*0.01,0,pk*10);
      m->SetLimitedVariable(1,"A52",0.1*pk,pk*0.001,0,pk*5);
      m->SetLimitedVariable(2,"B",50,1,0,1e7);
      m->SetFixedVariable(3,"a",av);
      m->SetLimitedVariable(4,"k",0.5,0.02,0.0,9.0);
      for(int g=0;g<4;g++) m->SetLimitedVariable(5+g,Form("A%d",g),20,0.5,0,1e5);
      m->SetLimitedVariable(npar-2,"mu",d.gsMu,0.005,d.gsMu-0.6,d.gsMu+0.6);
      m->SetLimitedVariable(npar-1,"sig",0.4,0.005,0.08,1.5);
      m->Minimize();
      if(m->MinValue()<bestL){ bestL=m->MinValue(); best.assign(m->X(),m->X()+npar); }
      delete m; }
    const double* p=best.data();
    printf("\n======== %s  theta_CM 20-25   a=%.2f k=%.2f  N=%.0f ========\n",
           CFG[ic].nm,p[3],p[4],d.N);
    printf("   %6s %8s %9s %8s   %8s %8s %8s  %s\n","Ex","data","model","pull","bkg","P1","P2","");
    for(double X=3.5;X<=9.01;X+=0.2){
      int ib=d.h->FindBin(X); double n=d.h->GetBinContent(ib);
      double m=B1::model(X,p); if(m<1e-9) m=1e-9;
      double Xt=X-p[npar-2];
      double bg=std::max(0.0,p[2])*B1::shp(Xt,p+3)*(B1::ACC?std::max(0.0,B1::ACC->Eval(Xt)):1.0);
      double a1=0,a2=0; const double* aa=p+5;
      for(int r=0;r<RK::NR;r++) if(RK::GRP[r]>=0&&d.in[r]){
        double v=std::max(0.0,aa[RK::GRP[r]])*((RK::RIDX[r]>=0)?RK::FIXRAT[RK::RIDX[r]]:1.0)*d.T[r].Eval(X);
        if(RK::GRP[r]==0) a1+=v; else if(RK::GRP[r]==1) a2+=v; }
      double pull=(n-m)/std::sqrt(std::max(m,1.0));
      char bar[40]; int k2=0; int np=(int)std::min(16.0,std::fabs(pull)*2.5);
      for(int j=0;j<np;j++) bar[k2++]= pull>0?'+':'-'; bar[k2]=0;
      printf("   %6.2f %8.0f %9.1f %8.2f   %8.1f %8.1f %8.1f  %s\n",X,n,m,pull,bg,a1,a2,bar); }
    RK::ONLYBIN=-1; CFG[0]=b0; CFG[1]=b1; }
  VZLO=0; RK::PAIRS=false; BB_QUIET=true; binbuild();
}
