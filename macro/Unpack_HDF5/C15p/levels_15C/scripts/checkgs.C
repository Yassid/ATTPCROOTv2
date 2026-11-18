#include "resfitD.C"
// Data vs model across the FULL range, forward bins, numbers not pictures.
void checkgs(){
  setCM(10,5,60);
  RK::exLo=3.5; RK::exHi=9.0; RK::ONETERM=true; RK::BKG=RK::HALON;
  RK::TIE12=false; RK::PAIRS=true; RK::NPERBIN=false;
  for(int i=0;i<RK::NR;i++) RK::USE[i]=(i<6);
  Cfg b0=CFG[0], b1=CFG[1];
  for(int ic=0;ic<1;ic++){
    CFG[0]=b0; CFG[1]=b1; CFG[ic].vzHi=900; CFG[ic].chi2=5.0;
    VZLO=300; BB_QUIET=true; binbuild();
    rkRun(CFG[ic],NCM,false);
    for(int b : {0,1,2}){
      RK::Bin& d=RK::B[b];
      printf("\n==== %s  theta_CM %.0f-%.0f  (gsAmp %.0f  gsMu %+.3f  FWHM %.2f/%.2f) ====\n",
             CFG[ic].nm,d.cmLo,d.cmHi,d.gsAmp,d.gsMu,2.355*d.sgL,2.355*d.sgR);
      printf("   %6s %9s %9s %8s   %8s %8s %8s %8s\n",
             "Ex","data","model","d/m","gs","5/2+","cont","P1+P2");
      for(double X=-1.0;X<=9.01;X+=0.5){
        int ib=d.h->FindBin(X); double n=d.h->GetBinContent(ib);
        double Xt=X-d.gsMu;
        double gs=d.gsAmp*shapeGT(X,d.gsMu,d.sgL,d.sgR,d.alL,d.alR);
        double e5=d.e5Amp*shapeGT(X,d.gsMu+0.740,d.sgL,d.sgR,d.alL,d.alR);
        double bg=RK::LASTBG[b]*RK::bshape(Xt,RK::LASTP0)*(d.ga?std::max(0.0,d.ga->Eval(Xt)):1.0);
        double rs=0; for(int r=0;r<RK::NR;r++) if(RK::GRP[r]>=0&&d.in[r]) rs+=RK::LASTA[b][r]*d.T[r].Eval(X);
        double m=gs+e5+bg+rs;
        printf("   %6.1f %9.0f %9.1f %8.2f   %8.1f %8.1f %8.1f %8.1f%s\n",
               X,n,m,(m>0?n/m:0),gs,e5,bg,rs, (X<RK::exLo)?"   [not fitted]":""); }
    }
  }
  CFG[0]=b0; CFG[1]=b1; VZLO=0; RK::PAIRS=false; BB_QUIET=true; binbuild();
}
