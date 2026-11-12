#include "binfit.C"
// Upper-limit scan.  Restricting the window keeps the fit out of the high-Ex region where the
// continuum model is poor and there is no known structure -- but a resonance must be DROPPED
// once the window no longer contains it, or its amplitude is set by a tail alone.
void exhi(){
  for(auto&k:CFG){
    printf("\n================ %s ================\n",k.nm);
    printf("  window        nBW   chi2/ndf per theta_CM bin 0..3        median\n");
    for(double hi : {7.0,8.0,9.0,10.0,12.0,14.0}){
      USE_LEBKG=false;
      REFTOT=0; for(int i=0;i<NB;i++) REF[i]=0;
      { double sa[NB]={0},tot=0;
        for(int ib=0;ib<6;ib++){ BR o=fitBin(k,ib,-2.0,hi,false,false);
          if(o.N<400) continue; for(int i=0;i<NB;i++) sa[i]+=o.A[i]; tot+=o.N; }
        for(int i=0;i<NB;i++) REF[i]=sa[i]; REFTOT=tot; }
      std::vector<double> cs;
      for(int ib=0;ib<4;ib++){ BR o=fitBin(k,ib,-2.0,hi,false,true);
        if(o.N>400) cs.push_back(o.chi2ndf); }
      int nin=NBW_IN;
      std::vector<double> sorted=cs; std::sort(sorted.begin(),sorted.end());
      printf("  [-2, %4.1f]     %d   ",hi,nin);
      for(double c:cs) printf(" %7.2f",c);
      printf("      %6.2f\n",sorted.empty()?-1:sorted[sorted.size()/2]);
    }
  }
}
