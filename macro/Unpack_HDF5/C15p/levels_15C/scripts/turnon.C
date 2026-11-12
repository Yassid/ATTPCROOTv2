#include "binfit.C"
// Scan the turn-on centre of the low-energy background.
//   minimum at Sn = 1.218         -> threshold physics (14C+n related)
//   minimum tracking the g.s. mu  -> part of the detector response
void turnon(){
  for(auto&k:CFG){
    REFTOT=0; for(int i=0;i<NB;i++) REF[i]=0;
    { double sa[NB]={0},tot=0; USE_LEBKG=false;
      for(int ib=0;ib<6;ib++){ BR o=fitBin(k,ib,-2.0,12.0,false,false);
        if(o.N<400) continue; for(int i=0;i<NB;i++) sa[i]+=o.A[i]; tot+=o.N; }
      for(int i=0;i<NB;i++) REF[i]=sa[i]; REFTOT=tot; }
    USE_LEBKG=true; FREE_TURNON=true;
    printf("\n== %s   -2lnL vs the background turn-on centre (bin 0; Sn = 1.218)\n",k.nm);
    printf("   centre   -2lnL   chi2    B      kappa |  BW significances (4.43 4.93 5.98 6.68 7.38 8.28)\n");
    double best=1e30,bc=0;
    for(double c=-1.0;c<=6.01;c+=0.5){ TURNON=c;
      BR o=fitBin(k,0,-2.0,12.0,false,true);
      printf("   %5.1f %8.0f %6.2f %7.1f %7.3f |",c,o.nll,o.chi2ndf,o.bgB,o.bgK);
      for(int i=0;i<NB;i++) printf(" %5.1f",o.err[i]>0?o.amp[i]/o.err[i]:0.0);
      printf("\n");
      if(o.nll<best){best=o.nll;bc=c;} }
    printf("   best turn-on %.2f MeV   (Sn = 1.218,  g.s. mu = %.3f)\n",bc,GSMU);
    FREE_TURNON=false; USE_LEBKG=false;
  }
}
