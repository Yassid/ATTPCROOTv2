#include "binfit.C"
// ADOPTED CONFIGURATION
//   Ex window [-2, 9]       keeps all six (p,p') resonances inside, and stops before the
//                           high-Ex region where the continuum model is poor and nothing is known
//   background turn-on 2.5  ABOVE the g.s., so its rising edge no longer sits under the g.s. tail
// With the turn-on at Sn the background had to choose between fitting the g.s. region and the
// resonance region; moved above the peak it fits both.
void final9(){
  for(auto&k:CFG){
    USE_LEBKG=true; FREE_TURNON=true; TURNON=2.5;
    REFTOT=0; for(int i=0;i<NB;i++) REF[i]=0;
    { double sa[NB]={0},tot=0; USE_LEBKG=false;
      for(int ib=0;ib<6;ib++){ BR o=fitBin(k,ib,-2.0,9.0,false,false);
        if(o.N<400) continue; for(int i=0;i<NB;i++) sa[i]+=o.A[i]; tot+=o.N; }
      for(int i=0;i<NB;i++) REF[i]=sa[i]; REFTOT=tot; USE_LEBKG=true; }
    printf("\n================ %s   Ex [-2, 9], background turn-on 2.5 ================\n",k.nm);
    printf(" bin  cm      N    chi2/ndf   B     kappa |");
    for(int i=0;i<NB;i++) printf("  BW%d(%.2f)",i+1,EPP[i]-PPGS);
    printf("\n");
    for(int ib=0;ib<6;ib++){
      BR o=fitBin(k,ib,-2.0,9.0,false,true);
      if(o.N<400){ printf(" %2d  %2.0f-%2.0f %6.0f   (too few)\n",ib,o.cmLo,o.cmHi,o.N); continue; }
      printf(" %2d  %2.0f-%2.0f %6.0f %8.2f %7.1f %6.3f |",ib,o.cmLo,o.cmHi,o.N,o.chi2ndf,o.bgB,o.bgK);
      for(int i=0;i<NB;i++){ double sg=o.err[i]>0?o.amp[i]/o.err[i]:0;
        printf(" %7.1f%s",sg,sg>5?"*":" "); }
      printf("\n"); }
    USE_LEBKG=false; FREE_TURNON=false;
  }
}
