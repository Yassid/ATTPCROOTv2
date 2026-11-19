#include "fitshow.C"
void look(){
  fitshow();
  for(int ic=0;ic<1;ic++) for(int b : {0,8,9}){
    Cfg b0=CFG[0], b1=CFG[1];
    CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; if(ic==0) CFG[ic].chi2=5.0;
    VZLO=300; BB_QUIET=true; binbuild();
    RK::SEED.clear(); rkRun(CFG[ic],NCM,false);
    for(int r=0;r<RK::NRAT;r++) RK::FIXRAT[r]=RK::LASTP[1+NCM*RK::nper()+r];
    RK::Bin& d=RK::B[b]; B1::BIN=b; B1::PEAKFREE=false; B1::MUFIX=d.gsMu;
    B1::NPAR=7; B1::NSH=2; B1::FORM=3; B1::EB=CFG[ic].eb;
    Pars A=parsOf(CFG[ic]); A.vzLo=300; Win w{14,86,CFG[ic].keLo,CFG[ic].keHi};
    B1::ACC=new TGraph();
    for(int i=0;i<=460;i++){ double X=i*0.05;
      B1::ACC->SetPoint(i,X,accBin(A,CFG[ic].eb,X,d.cmLo,d.cmHi,w)); }
    const double* p=QQ[ic][b];
    printf("\n==== %s  theta_CM %.0f-%.0f   N=%.0f   a=%.2f k=%.2f  chi2=%.2f ====\n",
           CFG[ic].nm,d.cmLo,d.cmHi,d.N,p[1],p[2],QX[ic][b]);
    printf("   in[] :"); for(int r=0;r<RK::NR;r++)
      printf(" %4.2f=%s",RK::EPP6[r]-RK::PPGS6,d.in[r]?"Y":"n"); printf("\n");
    printf("   s    :"); for(int g=0;g<4;g++) printf(" %8.3f",p[3+g]); printf("\n");
    printf("   amp  :"); for(int g=0;g<4;g++) printf(" %8.2f",p[3+g]*p[3+g]); printf("\n");
    printf("   %6s %8s %9s %8s | %8s %7s %7s %7s\n","Ex","data","model","pull","bkg","P1","P2","7.4/8.3");
    for(double X=2.0;X<=10.01;X+=0.5){
      int ib=d.h->FindBin(X); double n=d.h->GetBinContent(ib);
      double m=B1::model(X,p); if(m<1e-9) m=1e-9;
      double Xt=X-d.gsMu;
      double bg=std::max(0.0,p[0])*B1::shp(Xt,p+1)*(B1::ACC?std::max(0.0,B1::ACC->Eval(Xt)):1.0);
      double v[4]={0,0,0,0}; const double* a2=p+3;
      for(int r=0;r<RK::NR;r++) if(RK::GRP[r]>=0&&d.in[r])
        v[RK::GRP[r]]+=a2[RK::GRP[r]]*a2[RK::GRP[r]]
                       *((RK::RIDX[r]>=0)?RK::FIXRAT[RK::RIDX[r]]:1.0)*d.T[r].Eval(X);
      printf("   %6.1f %8.0f %9.1f %8.2f | %8.1f %7.1f %7.1f %7.1f\n",
             X,n,m,(n-m)/std::sqrt(std::max(m,1.0)),bg,v[0],v[1],v[2]+v[3]); }
    CFG[0]=b0; CFG[1]=b1; }
}
