#include "simfit.C"
// HOW MANY RESONANCES DOES THE DATA SUPPORT?
// Forward selection on the simultaneous fit: start with none, add the resonance that lowers
// -2lnL most, repeat.  Each addition costs 5 parameters (one amplitude per theta_CM bin), so a
// gain below ~11 (chi2 at 5 dof, p = 0.05) is not worth having.
void nres(){
  for(auto&k:CFG){
    printf("\n================ %s : forward selection ================\n",k.nm);
    for(int i=0;i<6;i++) SF::USEBW[i]=false;
    double prev=-1;
    printf("   step  added        -2lnL      gain   (5 par each; threshold ~11)\n");
    for(int step=1;step<=6;step++){
      int best=-1; double bestv=1e30;
      for(int i=0;i<6;i++){
        if(SF::USEBW[i]) continue;
        SF::USEBW[i]=true;
        double v=sfRun(k,false,false);
        if(v<bestv){bestv=v;best=i;}
        SF::USEBW[i]=false; }
      if(best<0) break;
      SF::USEBW[best]=true;
      double gain = (prev<0)?0:prev-bestv;
      printf("    %d    BW%d (%.2f)   %9.1f   %+8.1f%s\n",step,best+1,EPP[best]-PPGS,bestv,gain,
             (prev>0&&gain<11)?"   <- not worth it":"");
      prev=bestv; }
    for(int i=0;i<6;i++) SF::USEBW[i]=true;
  }
}
