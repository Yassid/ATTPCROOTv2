#include "binfit.C"
// HYPOTHESIS: the low-energy background is the low-energy TAIL OF THE ELASTIC PEAK (straggling
// and reconstruction loss at deuteron KE ~ 1.4 MeV), not a separate process.  Evidence so far:
// the excess events are deuterons with the same PID and track quality as the g.s., they sit only
// in the low-KE bin, and they turn on near threshold.
// SHARP PREDICTION: the background amplitude B must scale with the GROUND-STATE amplitude.
void bscale(){
  for(auto&k:CFG){
    REFTOT=0; for(int i=0;i<NB;i++) REF[i]=0;
    { double sa[NB]={0},tot=0; USE_LEBKG=false;
      for(int ib=0;ib<6;ib++){ BR o=fitBin(k,ib,-2.0,12.0,false,false);
        if(o.N<400) continue; for(int i=0;i<NB;i++) sa[i]+=o.A[i]; tot+=o.N; }
      for(int i=0;i<NB;i++) REF[i]=sa[i]; REFTOT=tot; }
    USE_LEBKG=true;
    printf("\n== %s   is B proportional to the g.s. amplitude?\n",k.nm);
    printf("   bin  cm      N      gs amp      B        B/gsAmp     w     kappa   <KE> of the bin\n");
    for(int ib=0;ib<5;ib++){
      BR o=fitBin(k,ib,-2.0,12.0,false,true);
      if(o.N<400) continue;
      // KE of this bin at Ex ~ 0 (two-body, exact)
      Pars A=parsOf(k); double thl,kk; labOf(A,k.eb,0.0,o.cmLo+TCMSTEP/2,thl,kk);
      printf("   %2d  %2.0f-%2.0f %7.0f %10.0f %8.1f %10.5f %7.2f %7.3f %10.2f MeV\n",
             ib,o.cmLo,o.cmHi,o.N,o.gsA,o.bgB,o.gsA>0?o.bgB/o.gsA:0,o.bgW,o.bgK,kk); }
    USE_LEBKG=false;
  }
}
