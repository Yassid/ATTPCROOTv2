#include "resfitB.C"
// Which functional form does the resonance-region background actually take?
// Compare by -2lnL at EQUAL parameter count (every form has exactly one shared shape parameter),
// and for the halo-breakup forms let the binding energy float: if it returns 15C's actual
// Sn = 1.218 MeV, the background is direct breakup and the shape is not empirical at all.
void bkform(){
  setCM(10,5,60);
  RG::exLo=3.0; RG::exHi=9.0; RG::TIE12=false;
  bool use[RG::NR]={true,true,true,true,true,true};
  for(int i=0;i<RG::NR;i++) RG::USE[i]=use[i];
  const char* nm[4]={"EXP   Bg*exp(-k*Ex)","SWAVE Bg*sqrt(Ex-Sb)/Ex^2","PWAVE Bg*(Ex-Sb)^1.5/Ex^3","LIN   Bg*(1-k*Ex)"};
  for(auto&k:CFG){
    printf("\n================ %s : background functional form ================\n",k.nm);
    printf("   %-28s %10s %9s   %s\n","form","-2lnL","shape par","");
    double best=1e18; int ib=-1; double v[4],sp[4];
    for(int m=0;m<4;m++){ RG::BKG=m; RG::SBFIX=false;
      v[m]=rgRun(k,NCM,false); sp[m]=RG::LASTP0;
      if(v[m]<best){best=v[m];ib=m;} }
    for(int m=0;m<4;m++)
      printf("   %-28s %10.1f %9.4f   %s\n",nm[m],v[m],sp[m],(m==ib)?"<- best":"");
    // s-wave with Sb pinned to theknown Sn
    RG::BKG=RG::SWAVE; RG::SBFIX=true; RG::SBVAL=1.218;
    double vfix=rgRun(k,NCM,false); RG::SBFIX=false;
    printf("   %-28s %10.1f %9.4f   (zero free shape parameters)\n",
           "SWAVE with Sb = Sn = 1.218",vfix,1.218);
    printf("   Delta(-2lnL) vs best = %+.1f for ONE fewer parameter\n",vfix-best);
  }
}
