#include "binfit.C"
// (A) nothing extra   (B) erf low-energy background   (C) a long tail on the elastic peak
void tailtest(){
  for(auto&k:CFG){
    REFTOT=0; for(int i=0;i<NB;i++) REF[i]=0;
    { double sa[NB]={0},tot=0; USE_LEBKG=false; USE_LONGTAIL=false;
      for(int ib=0;ib<6;ib++){ BR o=fitBin(k,ib,-2.0,12.0,false,false);
        if(o.N<400) continue; for(int i=0;i<NB;i++) sa[i]+=o.A[i]; tot+=o.N; }
      for(int i=0;i<NB;i++) REF[i]=sa[i]; REFTOT=tot; }
    printf("\n================ %s ================\n",k.nm);
    printf("  bin    N     (A) none        (B) erf bkg          (C) long elastic tail\n");
    printf("              chi2   -2lnL    chi2   d(-2lnL)     chi2   d(-2lnL)   frac   L(MeV)\n");
    for(int ib=0;ib<3;ib++){
      USE_LEBKG=false; USE_LONGTAIL=false; BR a=fitBin(k,ib,-2.0,12.0,false,true);
      if(a.N<400) continue;
      USE_LEBKG=true;  USE_LONGTAIL=false; BR b=fitBin(k,ib,-2.0,12.0,false,true);
      USE_LEBKG=false; USE_LONGTAIL=true;  BR c=fitBin(k,ib,-2.0,12.0,false,true);
      printf("  %2d %6.0f %7.2f %8.1f %7.2f %9.1f %9.2f %9.1f %7.4f %6.2f\n",
             ib,a.N,a.chi2ndf,a.nll,b.chi2ndf,a.nll-b.nll,c.chi2ndf,a.nll-c.nll,c.bgB,c.bgW); }
    USE_LEBKG=false; USE_LONGTAIL=false;
  }
}
