#include "resfit.C"
// Which resonances does (d,d') actually POPULATE?  A state seen in (p,p') need not be, and
// carrying one that is not just lets it absorb what the continuum misses.
// Forward selection in the Ex 3-9 window with 5 deg bins: each resonance costs 10 parameters
// (one amplitude per angular bin), so a gain below ~18 (chi2 at 10 dof, p=0.05) is not worth it.
void rsel(){
  setCM(10,5,60);
  for(auto&k:CFG){
    printf("\n================ %s : which resonances are populated ================\n",k.nm);
    for(int i=0;i<RF::NR;i++) RF::USE[i]=false;
    double prev=-1;
    printf("   step   added        -2lnL       gain    (10 par each; threshold ~18)\n");
    for(int step=1;step<=6;step++){
      int best=-1; double bestv=1e30;
      for(int i=0;i<RF::NR;i++){ if(RF::USE[i]) continue;
        RF::USE[i]=true; double v=rfRun(k,10,false);
        if(v<bestv){bestv=v;best=i;} RF::USE[i]=false; }
      if(best<0) break;
      RF::USE[best]=true;
      double gain=(prev<0)?0:prev-bestv;
      printf("    %d    %5.2f MeV   %9.1f   %+8.1f%s\n",step,RF::EPP6[best]-RF::PPGS6,bestv,gain,
             (prev>0&&gain<18)?"   <- STOP, not worth it":"");
      prev=bestv; }
    for(int i=0;i<RF::NR;i++) RF::USE[i]=true;
  }
}
